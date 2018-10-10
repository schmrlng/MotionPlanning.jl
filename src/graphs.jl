export NearNeighborGraph

# NeighborIterator
## Hacky, more performant alternative to Iterators.flatten((neighborhood, (D[v] for D in extras if v in keys(D))))
## Currently, to ensure iterate returns a Union{Nothing,XXX}, if isempty(neighborhood) then no extras will be considered
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
    iter_result === nothing ? nothing : (iter_result..., 0)
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

struct NearNeighborGraph{T,SS<:SampleSet,BVP<:SteeringBVP,NC<:NeighborhoodCache,DS<:NearNeighborDataStructure,D}
    nodes::SS
    bvp::BVP
    edge_cache::NC
    NN_data_structure::DS
    init_neighbors::Vector{D}
    goal_neighbors::Vector{D}
end
function NearNeighborGraph{T}(nodes::SS, bvp::BVP, edge_cache::NC, NN_data_structure::DS,
                              init_neighbors::Vector{D}, goal_neighbors::Vector{D}) where {T,SS,BVP,NC,DS,D}
    NearNeighborGraph{T,SS,BVP,NC,DS,D}(nodes, bvp, edge_cache, NN_data_structure, init_neighbors, goal_neighbors)
end
function NearNeighborGraph(nodes::SampleSet,
                           bvp::SteeringBVP,
                           edge_cache::NeighborhoodCache=MemoizedNC(nodes, bvp),
                           NN_data_structure::NearNeighborDataStructure=NullNNDS();
                           include_controls::BoolVal=Val(includes_controls(edge_cache)))
    T = neighbor_info_type(nodes, bvp, include_controls=include_controls)
    D = Dict{indextype(nodes),T}
    init_neighbors = Vector{D}(undef, 1)
    goal_neighbors = Vector{D}(undef, 0)
    NearNeighborGraph{T}(nodes, bvp, edge_cache, NN_data_structure, init_neighbors, goal_neighbors)
end
Base.getindex(G::NearNeighborGraph, i) = G.nodes[i]
# Base.getindex(G::NearNeighborGraph, i, j) = SteeringEdge
function setinit!(G::NearNeighborGraph{NeighborInfo{X,D,U}}, x::State, r) where {X,D,U}
    G.nodes.init = x    # TODO: also possibly reset!(G.edge_cache); add r as a field probably (RRT = nothing)
    G.init_neighbors[1] = reverse_neighbor_dict(X(0), neighbors(G, X(0), r, dir=Val(:F)), NeighborInfo{X,D,U})
end
function setgoal!(G::NearNeighborGraph{NeighborInfo{X,D,U}}, f::Base.Callable, r; count=1) where {X,D,U}
    resize!(G.nodes.goal_samples, count)
    resize!(G.goal_neighbors, count)
    for i in 1:count
        G.nodes.goal_samples[i] = f()
        G.goal_neighbors[i] = reverse_neighbor_dict(X(-i), neighbors(G, X(-i), r, dir=Val(:B)), NeighborInfo{X,D,U})
    end
end
setgoal!(G::NearNeighborGraph, x, r) = setgoal!(G, () -> x, r, count=1)
neighbor_info_type(G::NearNeighborGraph{T}) where {T} = T
indextype(G::NearNeighborGraph) = indextype(G.nodes)
# costtype(G::NearNeighborGraph) = costtype(neighbor_info_type(G))
# controlstype(G::NearNeighborGraph) = controlstype(neighbor_info_type(G))
includes_controls(G::NearNeighborGraph) = includes_controls(neighbor_info_type(G))

function neighbors(G::NearNeighborGraph{T}, v::Int, r; dir::F_or_B=Val(:F)) where {T}
    neighborhood = neighbors(G.edge_cache, G.NN_data_structure, G.nodes, G.bvp, v, r,
                             dir=dir, include_controls=Val(includes_controls(T)))
    NeighborIterator{T}(v, neighborhood, dir === Val(:F) ? G.goal_neighbors : G.init_neighbors)
    # if dir === Val(:F)
    #     Iterators.flatten((neighborhood, (nbhd[v] for nbhd in G.goal_neighbors if v in keys(nbhd))))
    # else
    #     Iterators.flatten((neighborhood, (nbhd[v] for nbhd in G.init_neighbors if v in keys(nbhd))))
    # end
end
function neighbors(G::NearNeighborGraph{T,<:TiledSampleSet}, v::TiledIndex, r; dir::F_or_B=Val(:F)) where {T}
    on, off = v.on_lattice, v.off_lattice
    vz = zero_lattice_component(v)
    lattice_neighborhood = neighbors(G.edge_cache, G.NN_data_structure, G.nodes, G.bvp, vz, r,
                                     dir=dir, include_controls=Val(includes_controls(T)))
    neighborhood = (shift_lattice_component(n, on) for n in lattice_neighborhood)
    NeighborIterator{T}(v, neighborhood, dir === Val(:F) ? G.goal_neighbors : G.init_neighbors)
    # if dir === Val(:F)
    #     Iterators.flatten((neighborhood, (nbhd[v] for nbhd in G.goal_neighbors if v in keys(nbhd))))
    # else
    #     Iterators.flatten((neighborhood, (nbhd[v] for nbhd in G.init_neighbors if v in keys(nbhd))))
    # end
end
function neighbors(f::Function, G::NearNeighborGraph, v::SampleIndex, r; dir::F_or_B=Val(:F))
    Iterators.filter(f, neighbors(G, v, r, dir=dir))
end
