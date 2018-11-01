export FLANNNNDS
const FLANNIndex = FLANN.FLANNIndex

mutable struct FLANNNNDS{T} <: MetricNNDS
    index::FLANNIndex{T}
    dataptr::Ptr{T}
    inds::Vector{Cint}
    dists::Vector{T}
end
function FLANNNNDS(nodes::ExplicitSampleSet{S}, bvp::GeometricSteering; max_neighbor_count=1000) where {S}
    T = eltype(S)
    params = FLANN.FLANNParameters(algorithm=FLANN.FLANN_INDEX_KDTREE_SINGLE, checks=-1, trees=1)
    inds  = fill(Cint(0), max_neighbor_count)
    dists = fill(T(0), max_neighbor_count)
    nnds = FLANNNNDS(FLANN.flann(vecs2mat(nodes.V), params, FLANN.Euclidean()),
                     convert(Ptr{T}, pointer(nodes.V)), inds, dists)
    finalizer(x -> FLANN.close(x.index), nnds)
    nnds
end

default_NN_data_structure(nodes::ExplicitSampleSet, bvp::GeometricSteering) = FLANNNNDS(nodes, bvp)

# addstates!(nnds::FLANNNNDS, x::State) = FLANN.addpoints!(nnds.index, collect_if_not_Array(x))    # does ref persist?
# addstates!(nnds::FLANNNNDS, xs::AbstractVector{<:State}) = FLANN.addpoints!(nnds.index, vecs2mat(xs))
function addstates!(nnds::FLANNNNDS{T}, v::SubArray{<:State})  where {T}   # assumes view v is at the end of its parent
    ptr = convert(Ptr{T}, pointer(v.parent))
    if ptr == nnds.dataptr
        FLANN.addpoints!(nnds.index, vecs2mat(v))
    else
        FLANN.close(nnds.index)
        params = FLANN.FLANNParameters(algorithm=FLANN.FLANN_INDEX_KDTREE_SINGLE, checks=-1, trees=1)
        nnds.index = FLANN.flann(vecs2mat(v.parent), params, FLANN.Euclidean())
        nnds.dataptr = ptr
    end
end

function neighbors(nnds::FLANNNNDS, nodes::ExplicitSampleSet, bvp::GeometricSteering, v::Int, r::Number;
                   dir::F_or_B=Val(:F), include_controls::BoolVal=Val(false))
    neighbors(nnds, nodes, bvp, nodes[v], r, dir=dir, include_controls=include_controls, omit=isequal(v))
end
function neighbors(nnds::FLANNNNDS, nodes::ExplicitSampleSet, bvp::GeometricSteering, x::State, r::Number;
                   dir::F_or_B=Val(:F), include_controls::BoolVal=Val(false), omit::F=always_false) where {F}
    near_indices = FLANN.inrange!(nnds.index, x, (bvp.constraints.b*r)^2, length(nnds.inds), nnds.inds, nnds.dists)[1]
    neighbor_info_iterator(nodes, bvp, x, Int.(near_indices), r, dir=dir, include_controls=include_controls, omit=omit)
end
function neighbors(nnds::FLANNNNDS, nodes::ExplicitSampleSet, bvp::GeometricSteering, v::Int, r::Nearest;
                   dir::F_or_B=Val(:F), include_controls::BoolVal=Val(false))
    neighbors(nnds, nodes, bvp, nodes[v], v > 0 ? Nearest(r.k + 1) : r,
              dir=dir, include_controls=include_controls, omit=isequal(v))
end
function neighbors(nnds::FLANNNNDS, nodes::ExplicitSampleSet, bvp::GeometricSteering, x::State, r::Nearest;
                   dir::F_or_B=Val(:F), include_controls::BoolVal=Val(false), omit::F=always_false) where {F}
    r.k > length(nnds.inds) && (resize!(nnds.inds, r.k); resize!(nnds.dists, r.k))
    near_indices = FLANN.knn!(nnds.index, x, r.k, nnds.inds, nnds.dists)[1]
    if r.k > length(nodes.V)    # == length(nnds.index)
        near_indices = view(near_indices, 1:length(nodes.V))
    elseif r.k < length(nnds.inds)
        near_indices = view(near_indices, 1:r.k)
    end
    neighbor_info_iterator(nodes, bvp, x, Int.(near_indices), dir=dir, include_controls=include_controls, omit=omit)
end
