export @any, @all

macro any(comp)
    if !isa(comp, Expr) || (comp.head != :comprehension)
        error("@any may only precede a comprehension")
    end
    compgen = comp.args[1]
    loop = Expr(:for, esc(compgen.args[2]),
                          quote
                              if $(esc(compgen.args[1]))
                                  result = true
                                  break
                              end
                          end)
    quote
        result = false
        $loop
        result
    end
end

macro all(comp)
    if !isa(comp, Expr) || (comp.head != :comprehension)
        error("@all may only precede a comprehension")
    end
    compgen = comp.args[1]
    loop = Expr(:for, esc(compgen.args[2]),
                          quote
                              if !$(esc(compgen.args[1]))
                                  result = false
                                  break
                              end
                          end)
    quote
        result = true
        $loop
        result
    end
end

@generated function blend{T1,T2}(b::StaticArray{Bool}, a1::StaticArray{T1}, a2::StaticArray{T2})
    if !(size(b) == size(a1) == size(a2))
        error("Dimensions must match. Got sizes $(size(b)), $(size(a1)), and $(size(a2))")
    end
    newtype = :(similar_type($a1, promote_type(T1, T2)))
    exprs = [:(b[$j] ? a1[$j] : a2[$j]) for j = 1:length(b)]
    return quote
        $(Expr(:meta, :inline))
        $(Expr(:call, newtype, Expr(:tuple, exprs...)))
    end
end

(::Type{SVector})(x::Vector) = SVector(x...)
(::Type{MVector})(x::Vector) = MVector(x...)

# Inspired by SparseVectors.jl
import Base: length, size, nnz, countnz
import Base.SparseArrays: nonzeros, nonzeroinds
export viewcol, nonzeroinds

typealias CVecView{T} SubArray{T,1,Vector{T},Tuple{UnitRange{Int64}},true}

immutable SparseVectorView{Tv,Ti<:Integer} <: AbstractSparseVector{Tv,Ti}
    n::Int                  # the number of elements
    nzind::CVecView{Ti}     # the indices of nonzeros
    nzval::CVecView{Tv}     # the values of nonzeros

    function SparseVectorView(n::Integer, nzind::CVecView{Ti}, nzval::CVecView{Tv})
        n >= 0 || throw(ArgumentError("The number of elements must be non-negative."))
        size(nzind) == size(nzval) || throw(DimensionMismatch("The lengths of nzind and nzval are inconsistent."))
        new(convert(Int, n), nzind, nzval)
    end
end
SparseVectorView{Tv,Ti}(n::Integer, nzind::CVecView{Ti}, nzval::CVecView{Tv}) = SparseVectorView{Tv,Ti}(n, nzind, nzval)

function viewcol(x::SparseMatrixCSC, j::Integer)
    1 <= j <= x.n || throw(BoundsError())
    r1 = convert(Int, x.colptr[j])
    r2 = convert(Int, x.colptr[j+1]) - 1
    rgn = r1:r2
    SparseVectorView(x.m, view(x.rowval, rgn), view(x.nzval, rgn))
end
length(x::SparseVectorView) = x.n
size(x::SparseVectorView) = (x.n,)
nnz(x::SparseVectorView) = length(x.nzval)
countnz(x::SparseVectorView) = countnz(x.nzval)
nonzeros(x::SparseVectorView) = x.nzval
nonzeroinds(x::SparseVectorView) = x.nzind


include("vec2Dutils.jl")