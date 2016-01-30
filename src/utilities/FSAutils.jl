import Base: hcat, full, pinv

function hcat{N,T}(X::Vec{N,T}...)
    result = Array(T, N, length(X))
    @inbounds for i in 1:length(X), j in 1:N
        result[j,i] = X[i][j]
    end
    result
end
full{N,T}(v::Vec{N,T}) = convert(Vector{T}, v)
full{M,N,T}(A::Mat{M,N,T}) = convert(Matrix{T}, A)
pinv{M,N,T}(A::Mat{M,N,T}) = Mat(pinv(full(A)))
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