export @any, @all

macro any(comp)
    if !isa(comp, Expr) || (comp.head != :comprehension)
        error("@any may only precede a comprehension")
    end
    loop = Expr(:for, esc(comp.args[2]),
                          quote
                              if $(esc(comp.args[1]))
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
    loop = Expr(:for, esc(comp.args[2]),
                          quote
                              if !$(esc(comp.args[1]))
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

macro unfix(f)
    if !isa(f, Expr) || (f.head != :function && f.head != :(=)) || !isa(f.args[1], Expr) || f.args[1].head != :call
        error("@unfix may only precede a function definition")
    end
    fLHS = copy(f.args[1])
    fRHS = copy(f.args[1])

    changed = false
    if isa(fRHS.args[1], Expr) && fRHS.args[1].head == :curly
        fRHS.args[1] = fRHS.args[1].args[1]
    end
    for i in 2:length(fLHS.args)
        farg = fLHS.args[i]
        if isa(farg, Expr) && farg.head == :(::)
            fRHS.args[i] = fRHS.args[i].args[1]
            argtype = farg.args[2]
            if argtype == :Vec
                farg.args[2] = :AbstractVector
                fRHS.args[i] = :(Vec($(fRHS.args[i])))
                changed = true
            elseif argtype == :Mat
                farg.args[2] = :AbstractMatrix
                fRHS.args[i] = :(Mat($(fRHS.args[i])))
                changed = true
            elseif isa(argtype, Expr) && argtype.head == :curly
                if argtype.args[1] == :Vec
                    if length(argtype.args) < 3
                        farg.args[2] = :AbstractVector
                    else
                        farg.args[2] = :(AbstractVector{$(argtype.args[3])})
                    end
                    fRHS.args[i] = :(Vec($(fRHS.args[i])))
                    changed = true
                elseif argtype.args[1] == :Mat
                    if length(argtype.args) < 4
                        farg.args[2] = :AbstractMatrix
                    else
                        farg.args[2] = :(AbstractMatrix{$(argtype.args[4])})
                    end
                    fRHS.args[i] = :(Mat($(fRHS.args[i])))
                    changed = true
                end
            end
        end
    end
    quote
        $(esc(f))
        if $changed
            $(esc(fLHS)) = $(esc(fRHS))
        end
    end
end

# Mat (just sticking this here; not sure where it really belongs)
dense{M,N,T}(A::Mat{M,N,T}) = convert(Matrix{T}, A)
changeprecision{M,N,T<:AbstractFloat}(::Type{T}, A::Mat{M,N}) = convert(Mat{M,N,T}, A)

# Inspired by SparseVectors.jl
import Base: length, size, nnz, countnz
import Base.SparseArrays: nonzeros, nonzeroinds
export subcol, nonzeroinds

typealias CVecSub{T} SubArray{T,1,Vector{T},Tuple{UnitRange{Int64}},true}

immutable SparseVectorSub{Tv,Ti<:Integer} <: AbstractSparseVector{Tv,Ti}
    n::Int                  # the number of elements
    nzind::CVecSub{Ti}     # the indices of nonzeros
    nzval::CVecSub{Tv}     # the values of nonzeros

    function SparseVectorSub(n::Integer, nzind::CVecSub{Ti}, nzval::CVecSub{Tv})
        n >= 0 || throw(ArgumentError("The number of elements must be non-negative."))
        nzind.dims == nzval.dims || throw(DimensionMismatch("The lengths of nzind and nzval are inconsistent."))
        new(convert(Int, n), nzind, nzval)
    end
end
SparseVectorSub{Tv,Ti}(n::Integer, nzind::CVecSub{Ti}, nzval::CVecSub{Tv}) = SparseVectorSub{Tv,Ti}(n, nzind, nzval)

function subcol(x::SparseMatrixCSC, j::Integer)
    1 <= j <= x.n || throw(BoundsError())
    r1 = convert(Int, x.colptr[j])
    r2 = convert(Int, x.colptr[j+1]) - 1
    rgn = r1:r2
    SparseVectorSub(x.m, sub(x.rowval, rgn), sub(x.nzval, rgn))
end
length(x::SparseVectorSub) = x.n
size(x::SparseVectorSub) = (x.n,)
nnz(x::SparseVectorSub) = length(x.nzval)
countnz(x::SparseVectorSub) = countnz(x.nzval)
nonzeros(x::SparseVectorSub) = x.nzval
nonzeroinds(x::SparseVectorSub) = x.nzind


include("vec2Dutils.jl")