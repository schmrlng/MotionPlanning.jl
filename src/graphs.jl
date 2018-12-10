export NearNeighborGraph
import DifferentialDynamicsModels: _dummy_iterate_state
import Base: RefValue

# Flatten2
## More performant alternative to, e.g., Iterators.flatten((neighborhood, (D[v] for D in extras if v in keys(D))))
struct Flatten2{T,I1,I2}
    a::I1
    b::I2
end
Flatten2{T}(a::I1, b::I2) where {T,I1,I2} = Flatten2{T,I1,I2}(a, b)
Base.IteratorSize(::Flatten2) = Base.SizeUnknown()
Base.eltype(::Type{<:Flatten2{T}}) where {T} = T
@inline function _start(itr::Flatten2)
    a_iter_result = iterate(itr.a)
    a_iter_result !== nothing && return (a_iter_result[1], 1, (a_iter_result[2], _dummy_iterate_state(itr.b)))
    b_iter_result = iterate(itr.b)
    b_iter_result !== nothing && return (b_iter_result[1], 2, (_dummy_iterate_state(itr.a), b_iter_result[2]))
    nothing
end
@inline function Base.iterate(itr::Flatten2, state=_start(itr))
    (state === nothing || state[2] > 2) && return nothing
    if state[2] == 1    # itr.a not done yet
        a_iter_result = iterate(itr.a, state[3][1])
        a_iter_result !== nothing && return (state[1], (a_iter_result[1], 1, (a_iter_result[2], state[3][2])))
        b_iter_result = iterate(itr.b)
    else
        b_iter_result = iterate(itr.b, state[3][2])
    end
    b_iter_result !== nothing && return (state[1], (b_iter_result[1], 2, (state[3][1], b_iter_result[2])))
    (state[1], (state[1], 3, state[3]))
end

struct NearNeighborGraph{T,SS<:SampleSet,BVP<:SteeringBVP,R,NC<:NeighborhoodCache,DS<:NearNeighborDataStructure,D}
    nodes::SS
    bvp::BVP
    r::R
    edge_cache::NC
    NN_data_structure::DS
    init_neighbors::Vector{D}
    goal_neighbors::Vector{D}
end
function NearNeighborGraph{T}(nodes::SS, bvp::BVP, r::R, edge_cache::NC, NN_data_structure::DS,
                              init_neighbors::Vector{D}, goal_neighbors::Vector{D}) where {T,SS,BVP,R,NC,DS,D}
    NearNeighborGraph{T,SS,BVP,R,NC,DS,D}(nodes, bvp, r, edge_cache, NN_data_structure, init_neighbors, goal_neighbors)
end
function NearNeighborGraph(nodes::SampleSet,
                           bvp::SteeringBVP;
                           near_style::Union{Val{:variable},Val{:fixedradius},Val{:fixedknn}}=Val(:fixedradius),
                           include_controls::BoolVal=Val(!(bvp isa GeometricSteering) && nodes isa TiledSampleSet),
                           edge_cache::NeighborhoodCache=default_edge_cache(nodes, bvp, include_controls, near_style),
                           NN_data_structure::NearNeighborDataStructure=default_NN_data_structure(nodes, bvp,
                                                                                                  include_controls))
    T = neighbor_info_type(nodes, bvp, include_controls=include_controls)
    X = indextype(T)
    R = costtype(T)
    D = Dict{X,T}
    r = near_style === Val(:variable)    ? nothing :
        near_style === Val(:fixedradius) ? RefValue{R}(zero(R)) :
     #= near_style === Val(:fixedknn)   =# RefValue{Nearest}(Nearest(0))
    init_neighbors = D[D()]
    goal_neighbors = D[]
    NearNeighborGraph{T}(nodes, bvp, r, edge_cache, NN_data_structure, init_neighbors, goal_neighbors)
end
neighbor_info_type(G::NearNeighborGraph{T}) where {T} = T
for f in (:indextype, :costtype, :controlstype, :includes_controls)
    @eval $f(G::NearNeighborGraph{T}) where {T} = $f(T)
end

default_edge_cache(nodes, bvp, include_controls, ::Val{:variable}) = StreamingNC()
default_edge_cache(nodes, bvp, include_controls, ::Any) = MemoizedNC(nodes, bvp, include_controls=include_controls)
const default_NN_data_structure_type = Ref(:NullNNDS)    # probably not the most performant way to do this
function default_NN_data_structure(nodes, bvp, include_controls)
    getfield(@__MODULE__, default_NN_data_structure_type[])(nodes, bvp, include_controls=include_controls)
end

Base.length(G::NearNeighborGraph) = length(G.nodes)
Base.getindex(G::NearNeighborGraph, i) = G.nodes[i]
Base.getindex(G::NearNeighborGraph, i, j) = SteeringEdge(G.bvp, G[i], G[j])    # TODO: use cached controls if applicable

setradius!(G::NearNeighborGraph{<:Any,<:Any,<:Any,<:RefValue}, r) = (G.r[] != r && reset!(G.edge_cache); G.r[] = r)
setradius!(G::NearNeighborGraph{<:Any,<:Any,<:Any,Nothing},    r) = error("Cannot setradius! for variable radius graph")
function setinit!(G::NearNeighborGraph{NeighborInfo{X,D,U}}, x::State) where {X,D,U}
    G.nodes.init = x
    G.init_neighbors[1] = Dict{X,NeighborInfo{X,D,U}}()
end
function setgoal!(G::NearNeighborGraph{NeighborInfo{X,D,U}}, f; sample_count=1) where {X,D,U}
    resize!(G.nodes.goal_samples, sample_count)
    resize!(G.goal_neighbors, sample_count)
    for i in 1:sample_count
        G.nodes.goal_samples[i] = f()
        G.goal_neighbors[i] = Dict{X,NeighborInfo{X,D,U}}()
    end
end
setgoal!(G::NearNeighborGraph, x::State) = setgoal!(G, () -> x, sample_count=1)
function compute_init_and_goal_neighbors!(G::NearNeighborGraph{NeighborInfo{X,D,U}}) where {X,D,U}
    foreach(empty!, G.init_neighbors)
    foreach(empty!, G.goal_neighbors)
    G.init_neighbors[1] = reverse_neighbor_dict(X(0), neighbors(G, X(0), dir=Val(:F)), NeighborInfo{X,D,U})
    for i in eachindex(G.goal_neighbors)
        bvp_result = keep_controls_if(G.bvp(G[X(0)], G[X(-i)]), Val(includes_controls(G)))
        bvp_result.cost <= G.r[] && (G.init_neighbors[1][X(-i)] = (index=X(0), bvp_result...))
        G.goal_neighbors[i] = reverse_neighbor_dict(X(-i), neighbors(G, X(-i), dir=Val(:B)), NeighborInfo{X,D,U})
    end
end
function addstates!(G::NearNeighborGraph, X)
    # foreach(A -> addstates!(A, X), (G.nodes, G.edge_cache, G.NN_data_structure))
    Xview = addstates!(G.nodes, X)
    addstates!(G.edge_cache, Xview)
    addstates!(G.NN_data_structure, Xview)
end

function neighbors(G::NearNeighborGraph{T,<:ExplicitSampleSet}, v_or_x; dir::F_or_B=Val(:F)) where {T}
    neighborhood = neighbors(G.edge_cache, G.NN_data_structure, G.nodes, G.bvp, v_or_x, G.r[],
                             dir=dir, include_controls=Val(includes_controls(T)))
    include_init_and_goal_neighbors(G, neighborhood, v_or_x, dir)
end
function neighbors(G::NearNeighborGraph{T,<:ExplicitSampleSet}, v_or_x, r; dir::F_or_B=Val(:F)) where {T}
    neighborhood = neighbors(StreamingNC(), G.NN_data_structure, G.nodes, G.bvp, v_or_x, r,
                             dir=dir, include_controls=Val(includes_controls(T)))
    include_init_and_goal_neighbors(G, neighborhood, v_or_x, r, dir)
end
function neighbors(G::NearNeighborGraph{T,<:TiledSampleSet}, v::TiledIndex; dir::F_or_B=Val(:F)) where {T}
    relative_nbhd = neighbors(G.edge_cache, G.NN_data_structure, G.nodes, G.bvp, zero_lattice_component(v), G.r[],
                              dir=dir, include_controls=Val(includes_controls(T)))
    neighborhood = (shift_lattice_component(n, v.on_lattice) for n in relative_nbhd)
    include_init_and_goal_neighbors(G, neighborhood, v, dir)
end
function neighbors(f::F, G::NearNeighborGraph, v::SampleIndex; dir::F_or_B=Val(:F)) where {F}
    Iterators.filter(f, neighbors(G, v, dir=dir))
end
function neighbors(f::F, G::NearNeighborGraph, v::SampleIndex, r; dir::F_or_B=Val(:F)) where {F}
    Iterators.filter(f, neighbors(G, v, r, dir=dir))
end

function include_init_and_goal_neighbors(G::NearNeighborGraph, neighborhood, v_or_x, dir)
    include_init_and_goal_neighbors(G, neighborhood, v_or_x, G.r[], dir)
end
function include_init_and_goal_neighbors(G::NearNeighborGraph{T,<:SampleSet,<:SteeringBVP,<:RefValue{<:Number}},
                                         neighborhood, v::SampleIndex, dir) where {T}
    Flatten2{T}(neighborhood, (D[v] for D in (dir === Val(:F) ? G.goal_neighbors : G.init_neighbors) if v in keys(D)))
end
function include_init_and_goal_neighbors(G::NearNeighborGraph, neighborhood, v::SampleIndex, r, dir)
    include_init_and_goal_neighbors(G, neighborhood, G[v], r, dir)
end
function include_init_and_goal_neighbors(G::NearNeighborGraph{T}, neighborhood, x::State, r::Number, dir) where {T}
    extra_inds = dir === Val(:F) ? -eachindex(G.nodes.goal_samples) : 0:0    # actually doubt type stability matters
    extras = neighbor_info_iterator(G.nodes, G.bvp, x, extra_inds, r,
                                    dir=dir, include_controls=Val(includes_controls(T)))
    Flatten2{T}(neighborhood, Iterators.filter(n -> n.cost <= r, extras))
end
function include_init_and_goal_neighbors(G::NearNeighborGraph{T}, neighborhood, x::State, r::Nearest, dir) where {T}
    extra_inds = dir === Val(:F) ? -eachindex(G.nodes.goal_samples) : 0:0    # actually doubt type stability matters
    neighborhood = collect_if_not_Array(neighborhood)
    extras = neighbor_info_iterator(G.nodes, G.bvp, x, extra_inds, neighborhood[1].cost,
                                    dir=dir, include_controls=Val(includes_controls(T)))
    foreach(n -> update_knn!(neighborhood, n), extras)
    neighborhood
end

# Generic Tree Node Info (RRT, RRTstar)
const TreeNodeInfo{X,D} = NamedTuple{(:parent, :cost_to_come),Tuple{X,D}}
TreeNodeInfo{X,D}(; parent=zero(X), cost_to_come=zero(D)) where {X,D} = (parent=parent, cost_to_come=cost_to_come)

const FullTreeNodeInfo{X,D} = NamedTuple{(:parent, :cost_to_come, :children),Tuple{X,D,Vector{X}}}    # better name?
function FullTreeNodeInfo{X,D}(; parent=zero(X), cost_to_come=zero(D), children=Vector{X}(undef, 0)) where {X,D}
    (parent=parent, cost_to_come=cost_to_come, children=children)
end

# TODOs
# - add `dirty` flag to track NearNeighborGraph init_neighbors/goal_neighbors
# - also, enforce compute_init_and_goal_neighbors! if R<:RefValue{<:Number} -- otherwise, brittleness
# ✓ NeighborInfo U should maybe be <: Union{...,Missing} instead of Union{...,Nothing}
