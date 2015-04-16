import SymPy     # using "Using" gives me a weird diff conflict
using FastAnonymous

### Symbolic stuff

function expAt(A::Matrix, t::SymPy.Sym)
    n = size(A,1)
    maximum(abs(A^n)) > 0 && error("TODO: implement more cases than nilpotent A! (e.g. diagonalizable)")
    sum([A^i*t^i/factorial(i) for i in 0:n-1])
end

# Need this utility function since convert produces f : Sym -> Sym
function Sym2Function(s::SymPy.Sym, args::Union(Symbol, Expr) = :t, replace_rule = ())
    sliteral = string(SymPy.simplify(SymPy.expand(s)))    # sure, why not?
    length(replace_rule) > 0 && (sliteral = replace(sliteral, replace_rule...))
    @eval @anon $args -> $(parse(sliteral))
end
function Sym2Function(s::Vector{SymPy.Sym}, args::Union(Symbol, Expr) = :t, replace_rule = ())
    vliteral = "[" * join(map(x -> string(SymPy.simplify(SymPy.expand(x))), s), "; ") * "]"
    length(replace_rule) > 0 && (vliteral = replace(vliteral, replace_rule...))
    @eval @anon $args -> $(parse(vliteral))
end
function Sym2Function(s::Matrix{SymPy.Sym}, args::Union(Symbol, Expr) = :t, replace_rule = ())
    mliteral = "[" * join(mapslices(rw -> join(rw, " "), map(x -> string(SymPy.simplify(SymPy.expand(x))), s), 2), "; ") * "]"
    length(replace_rule) > 0 && (mliteral = replace(mliteral, replace_rule...))
    @eval @anon $args -> $(parse(mliteral))
end

### Time-optimal 2BVP solution

immutable LinearQuadratic2BVP
    Ginv::DataType      # weird side effect of FastAnonymous
    cost::DataType
    dcost::DataType
    ddcost::DataType
    x::DataType
end
function LinearQuadratic2BVP(A::Matrix, B::Matrix, c::Vector, R::Matrix)
    n = size(A, 1)
    t = SymPy.symbols("t", real=true)
    s = SymPy.symbols("s", real=true)
    xS = [SymPy.symbols(join(["x$i" for i in 1:n], ", "), real = true)...]
    yS = [SymPy.symbols(join(["y$i" for i in 1:n], ", "), real = true)...]
    replace_rule = (r"[xy]\d+", s -> s[1:1] * "[" * s[2:end] * "]")

    expAtS = expAt(A, t)
    expAsS = expAt(A, s)
    GS = SymPy.integrate(expAtS*B*inv(R)*B'*expAtS', t)
    GinvS = inv(GS)
    cdriftS = SymPy.integrate(expAtS, t)*c
    xbarS = expAtS*xS + cdriftS
    costS = t + ((yS - xbarS)'*GinvS*(yS - xbarS))[1]
    dcostS = SymPy.diff(costS, t)
    ddcostS = SymPy.diff(costS, t, 2)
    xS = expAsS*xS + SymPy.integrate(expAsS, s)*c +
         SymPy.integrate(expAsS*B*inv(R)*B'*expAsS', s)*expAt(A, t - s)'*GinvS*(yS-xbarS)
    
    LinearQuadratic2BVP(Sym2Function(GinvS, :(t), replace_rule),
                        Sym2Function(costS, :(x, y, t), replace_rule),
                        Sym2Function(dcostS, :(x, y, t), replace_rule),
                        Sym2Function(ddcostS, :(x, y, t), replace_rule),
                        Sym2Function(xS, :(x, y, t, s), replace_rule))
end

function topt_bisection{dc}(::Type{dc}, x0, x1, tm; tol = 1e-6)
    # Bisection
    b = tm
    dc(x0, x1, b) < 0 && return tm
    a = .01
    while dc(x0, x1, a) > 0; a /= 2; end
    m = 0.
    cdval = 1.
    while abs(cdval) > tol
        m = (a+b)/2
        cdval = dc(x0, x1, m)
        cdval > 0 ? b = m : a = m
    end
    m
end

function topt_newton{dc, ddc}(::Type{dc}, ::Type{ddc}, x0, x1, tm; tol = 1e-6)
    # Bisection / Newton's method combo
    b = tm
    dc(x0, x1, b) < 0 && return tm
    a = .01*tm
    while dc(x0, x1, a) > 0; a /= 2; end
    t = tm / 2
    cdval = dc(x0, x1, t)
    while abs(cdval) > tol
        t = t - cdval / ddc(x0, x1, t)
        (t < a || t > b) && (t = (a+b)/2)
        cdval = dc(x0, x1, t)
        cdval > 0 ? b = t : a = t
    end
    t
end

function steer(L::LinearQuadratic2BVP, x0, x1, r)
    t = topt_newton(L.dcost, L.ddcost, x0, x1, r)
    L.cost(x0, x1, t), t
end