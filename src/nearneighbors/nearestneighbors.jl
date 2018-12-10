export NNTreeNNDS, use_NNTreeNNDS, use_NearestNeighbors
const NNTree = NearestNeighbors.NNTree
const KDTree = NearestNeighbors.KDTree
const Euclidean = NearestNeighbors.Euclidean

struct NNTreeNNDS{T<:NNTree} <: MetricNNDS
    tree::T
end
NNTreeNNDS(nodes::ExplicitSampleSet, bvp::GeometricSteering; kwargs...) = NNTreeNNDS(KDTree(nodes.V, Euclidean()))

use_NNTreeNNDS() = (default_NN_data_structure_type[] = :NNTreeNNDS)
use_NearestNeighbors() = use_NNTreeNNDS()
use_NearestNeighbors()    # set to default upon `import NearestNeighbors`

function neighbors(nnds::NNTreeNNDS, nodes::ExplicitSampleSet, bvp::GeometricSteering, v::Int, r;
                   dir::F_or_B=Val(:F), include_controls::BoolVal=Val(false))
    neighbors(nnds, nodes, bvp, nodes[v], r, dir=dir, include_controls=include_controls, omit=isequal(v))
end
function neighbors(nnds::NNTreeNNDS, nodes::ExplicitSampleSet, bvp::GeometricSteering, x::State, r;
                   dir::F_or_B=Val(:F), include_controls::BoolVal=Val(false), omit::F=always_false) where {F}
    near_indices = NearestNeighbors.inrange(nnds.tree, x, bvp.constraints.b*r)
    neighbor_info_iterator(nodes, bvp, x, near_indices, r, dir=dir, include_controls=include_controls, omit=omit)
    # gen = dir === Val(:F) ? ((index=i, bvp(x, nodes[i], r)...) for i in near_indices if !omit(i)) :
    #                         ((index=i, bvp(nodes[i], x, r)...) for i in near_indices if !omit(i))
    # (keep_controls_if(n, include_controls) for n in gen)
end
function neighbors(nnds::NNTreeNNDS, nodes::ExplicitSampleSet, bvp::GeometricSteering, x::State, r::Nearest;
                   dir::F_or_B=Val(:F), include_controls::BoolVal=Val(false), omit::F=always_false) where {F}
    near_indices = NearestNeighbors.knn(nnds.tree, x, r.k, false, omit)[1]
    neighbor_info_iterator(nodes, bvp, x, near_indices, dir=dir, include_controls=include_controls)
end
