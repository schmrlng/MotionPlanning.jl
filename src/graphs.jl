export NearNeighborGraph
import DifferentialDynamicsModels: _dummy_iterate_state
import Base: RefValue

# NeighborIterator
## Hacky, more performant alternative to Iterators.flatten((neighborhood, (D[v] for D in extras if v in keys(D))))
struct NeighborIterator{T,X,N,D}
    v::X
    neighborhood::N
    extras::Vector{D}
end
function NeighborIterator{T}(v::X, neighborhood::N, extras::Vector{D}) where {T,X,N,D}
    NeighborIterator{T,X,N,D}(v, neighborhood, extras)
end
Base.IteratorSize(::NeighborIterator) = Base.SizeUnknown()
Base.eltype(::Type{<:NeighborIterator{T}}) where {T} = T
@inline function _start(itr::NeighborIterator)
    iter_result = iterate(itr.neighborhood)
    iter_result !== nothing && return (iter_result..., 0)
    for i in 1:length(itr.extras)
        itr.v in keys(itr.extras[i]) && return (itr.extras[i][itr.v], _dummy_iterate_state(itr.neighborhood), 1)
    end
    nothing
end
@inline function Base.iterate(itr::NeighborIterator, state=_start(itr))
    state === nothing && return nothing
    if state[3] == 0
        iter_result = iterate(itr.neighborhood, state[2])
        iter_result === nothing ? (state[1], (state[1], state[2], 1)) : (state[1], (iter_result..., 0))
    else
        for i in state[3]:length(itr.extras)
            itr.v in keys(itr.extras[i]) && return (itr.extras[i][itr.v], (state[1], state[2], i + 1))
        end
        nothing
    end
end
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
                           include_controls::BoolVal=Val(false),
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
indextype(G::NearNeighborGraph) = indextype(G.nodes)
costtype(G::NearNeighborGraph) = costtype(neighbor_info_type(G))
controlstype(G::NearNeighborGraph) = controlstype(neighbor_info_type(G))
includes_controls(G::NearNeighborGraph) = includes_controls(neighbor_info_type(G))

default_edge_cache(nodes, bvp, include_controls, ::Val{:variable}) = StreamingNC()
default_edge_cache(nodes, bvp, include_controls, ::Any) = MemoizedNC(nodes, bvp, include_controls=include_controls)
default_NN_data_structure(nodes, bvp, include_controls) = default_NN_data_structure(nodes, bvp)
default_NN_data_structure(nodes, bvp) = NullNNDS()

Base.getindex(G::NearNeighborGraph, i) = G.nodes[i]
# Base.getindex(G::NearNeighborGraph, i, j) = SteeringEdge

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
    r = G.r[]
    G.init_neighbors[1] = reverse_neighbor_dict(X(0), neighbors(G, X(0), r, dir=Val(:F)), NeighborInfo{X,D,U})
    for i in eachindex(G.goal_neighbors)
        bvp_result = keep_controls_if(G.bvp(G[X(0)], G[X(-i)]), Val(U !== Nothing))
        bvp_result.cost <= r && (G.init_neighbors[1][X(-i)] = (index=X(0), bvp_result...))
        G.goal_neighbors[i] = reverse_neighbor_dict(X(-i), neighbors(G, X(-i), r, dir=Val(:B)), NeighborInfo{X,D,U})
    end
end
addstates!(G::NearNeighborGraph, X) = nothing    # TODO

# function include_init_and_goal_neighbors()

function neighbors(G::NearNeighborGraph{T,<:ExplicitSampleSet,<:Any,<:RefValue}, v::Int;
                   dir::F_or_B=Val(:F)) where {T}
    neighborhood = neighbors(G.edge_cache, G.NN_data_structure, G.nodes, G.bvp, v, G.r[],
                             dir=dir, include_controls=Val(includes_controls(T)))
    # NeighborIterator{T}(v, neighborhood, dir === Val(:F) ? G.goal_neighbors : G.init_neighbors)
    Flatten2{T}(neighborhood, (D[v] for D in (dir === Val(:F) ? G.goal_neighbors : G.init_neighbors) if v in keys(D)))
end
function neighbors(G::NearNeighborGraph{T,<:ExplicitSampleSet,<:Any,RefValue{Nearest}}, v::Int;
                   dir::F_or_B=Val(:F)) where {T}
    neighborhood = neighbors(G.edge_cache, G.NN_data_structure, G.nodes, G.bvp, v, G.r[],
                             dir=dir, include_controls=Val(includes_controls(T)))
    extra_inds = dir === Val(:F) ? -eachindex(G.nodes.goal_samples) : 0:0    # actually doubt type stability matters
    extras = neighbor_info_iterator(G.nodes, G.bvp, G[v], extra_inds, neighborhood[1].cost,
                                    dir=dir, include_controls=Val(includes_controls(T)), omit=isequal(v))
    foreach(n -> update_knn!(neighborhood, n), extras)
    neighborhood
end
function neighbors(G::NearNeighborGraph{T,<:ExplicitSampleSet}, v::Int, r; dir::F_or_B=Val(:F)) where {T}
    neighborhood = neighbors(StreamingNC(), G.NN_data_structure, G.nodes, G.bvp, v, r,
                             dir=dir, include_controls=Val(includes_controls(T)))
    extra_inds = dir === Val(:F) ? -eachindex(G.nodes.goal_samples) : 0:0
    extras = neighbor_info_iterator(G.nodes, G.bvp, G[v], extra_inds, r isa Nearest ? neighborhood[1].cost : r,
                                    dir=dir, include_controls=Val(includes_controls(T)), omit=isequal(v))
    if r isa Nearest
        foreach(n -> update_knn!(neighborhood, n), extras)
        neighborhood
    else
        Flatten2{T}(neighborhood, Iterators.filter(n -> n.cost <= r, extras))
    end
end
# function neighbors(G::NearNeighborGraph{T,<:ExplicitSampleSet,<:Any,Nothing}, v::Int, r;
#                    dir::F_or_B=Val(:F)) where {T}
#     neighborhood = neighbors(G.edge_cache, G.NN_data_structure, G.nodes, G.bvp, v, r,
#                              dir=dir, include_controls=Val(includes_controls(T)))
#     extras = dir === Val(:F) ? ((index=-i, G.bvp(G[v], G[-i])...) for i in eachindex(G.nodes.goal_samples)) :
#                                ((index=i, G.bvp(G[i], G[v])...) for i in 0:0)    # actually doubt type stability matters
#     Flatten2{T}(neighborhood, (keep_controls_if(n, Val(includes_controls(T))) for n in extras if n.cost <= r))
# end
function neighbors(G::NearNeighborGraph{T,<:TiledSampleSet}, v::TiledIndex, r=G.r[]; dir::F_or_B=Val(:F)) where {T}
    on, off = v.on_lattice, v.off_lattice
    vz = zero_lattice_component(v)
    lattice_neighborhood = neighbors(G.edge_cache, G.NN_data_structure, G.nodes, G.bvp, vz, r,
                                     dir=dir, include_controls=Val(includes_controls(T)))
    neighborhood = (shift_lattice_component(n, on) for n in lattice_neighborhood)
    # NeighborIterator{T}(v, neighborhood, dir === Val(:F) ? G.goal_neighbors : G.init_neighbors)
    Flatten2{T}(neighborhood, (D[v] for D in (dir === Val(:F) ? G.goal_neighbors : G.init_neighbors) if v in keys(D)))
end
function neighbors(f::F, G::NearNeighborGraph, v::SampleIndex; dir::F_or_B=Val(:F)) where {F}
    Iterators.filter(f, neighbors(G, v, dir=dir))
end
function neighbors(f::F, G::NearNeighborGraph, v::SampleIndex, r; dir::F_or_B=Val(:F)) where {F}
    Iterators.filter(f, neighbors(G, v, r, dir=dir))
end

# Generic Tree Node Info (RRT, RRTstar)
const TreeNodeInfo{X,D} = NamedTuple{(:parent, :cost_to_come),Tuple{X,D}}
TreeNodeInfo{X,D}(; parent=zero(X), cost_to_come=zero(D)) where {X,D} = (parent=parent, cost_to_come=cost_to_come)
