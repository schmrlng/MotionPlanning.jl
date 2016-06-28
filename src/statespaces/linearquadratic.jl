export LinearQuadratic, LinearQuadratic2BVP
export LinearQuadraticQuasiMetricSpace, DoubleIntegrator
export steer, steer_pairwise

### QuasiMetric Typedef/Evaluation
immutable LinearQuadratic2BVP{T<:AbstractFloat,
                              FGinv<:Function,
                              FexpAt<:Function,
                              Fcdrift<:Function,
                              Fcost<:Function,
                              Fdcost<:Function,
                              Fddcost<:Function,
                              Fx<:Function}
    A::Matrix{T}
    B::Matrix{T}
    c::Vector{T}
    R::Matrix{T}
    Ginv::FGinv
    expAt::FexpAt
    cdrift::Fcdrift
    cost::Fcost
    dcost::Fdcost
    ddcost::Fddcost
    x::Fx
end
type LinearQuadratic{T<:AbstractFloat,L<:LinearQuadratic2BVP} <: QuasiMetric
    bvp::L
    cmax::T
end
LinearQuadratic{T}(A::Matrix{T}, B::Matrix{T}, c::Vector{T}, R::Matrix{T}, cmax::T = T(1)) =
    LinearQuadratic(LinearQuadratic2BVP(A, B, c, R), cmax)
setup_steering(d::LinearQuadratic, r) = (d.cmax = r)
steer(d::LinearQuadratic, v::Vec, w::Vec) = steer(d.bvp, v, w, d.cmax)

evaluate(d::LinearQuadratic, v::Vec, w::Vec) = steer(d, v, w)[1]
changeprecision{T<:AbstractFloat}(::Type{T}, d::LinearQuadratic) =
    LinearQuadratic(LinearQuadratic2BVP((map(T,x) for x in (d.bvp.A, d.bvp.B, d.bvp.c, d.bvp.R))...), T(d.cmax))

### QuasiMetric Space Instantiation
LinearQuadraticQuasiMetricSpace(lo::Vec, hi::Vec, A::Matrix, B::Matrix, c::Vector, R::Matrix, C::Matrix) =
    RealVectorStateSpace(lo, hi, LinearQuadratic(A, B, c, R), OutputMatrix(C))
function DoubleIntegrator(d::Int, lo = zeros(d), hi = ones(d); vmax = 1.5, r = 1.)
    A = [zeros(d,d) eye(d); zeros(d,2d)]
    B = [zeros(d,d); eye(d)]
    c = zeros(2d)
    R = r*eye(d)
    C = [eye(d) zeros(d,d)]
    LinearQuadraticQuasiMetricSpace(Vec([lo; -vmax*ones(d)]), Vec([hi; vmax*ones(d)]), A, B, c, R, C)
end
function WebbJvdB13quad10d()
    g = 9.8
    m = .5
    # l = 
    A32 = [0 g; -g 0; 0 0]
    A = [zeros(3,3) eye(3)      zeros(3,2)  zeros(3,2);
         zeros(3,3) zeros(3,3)  A32         zeros(3,2);
         zeros(2,3) zeros(2,3)  zeros(2,2)  eye(2);
         zeros(2,10)]
    B = [zeros(5, 3);
         ]
    c = zeros(10)
end

function helper_data_structures{S<:Vec}(V::Vector{S}, M::LinearQuadratic, batchsize = 1001)
    N = length(V)
    DUpairs = [steer_pairwise(M.bvp, V, V[rng], M.cmax) for rng in [i:min(i+batchsize-1, N) for i in 1:batchsize:N]]
    Dmat = hcat([D for (D,U) in DUpairs]...)
    Umat = hcat([U for (D,U) in DUpairs]...)
    DSF = BruteDistanceDS(Dmat')
    DSB = BruteDistanceDS(Dmat)
    US = EmptyControlCache()   # TODO: non-empty ControlDS
    DSF, US, DSB, US
end

### Steering
propagate(d::LinearQuadratic, v::Vec, u::DurationAndTargetControl) = u.x
propagate(d::LinearQuadratic, v::Vec, u::DurationAndTargetControl, s::AbstractFloat) =
    s <= 0 ? v : s >= u.t ? u.x : d.bvp.x(v, u.x, u.t, s)

steering_control(d::LinearQuadratic, v::Vec, w::Vec) = DurationAndTargetControl(steer(d,v,w)[2], w)
function collision_waypoints(d::LinearQuadratic, v::Vec, w::Vec)
    t = steer(d, v, w)[2]
    [d.bvp.x(v, w, t, s) for s in linspace(typeof(t)(0), t, 5)]
end

### LQ2BVP Nuts and Bolts
import SymPy

## Symbolic stuff
function expAt(A::Matrix, t::SymPy.Sym)
    n = size(A,1)
    maximum(abs(A^n)) > 0 && error("TODO: implement more cases than nilpotent A! (e.g. diagonalizable)")
    sum([A^i*(t^i/factorial(i)) for i in 0:n-1])
end
function Sym2Function(s::SymPy.Sym, args::Union{Symbol, Expr} = :t, replace_rules = ())
    sliteral = string(SymPy.simplify(SymPy.expand(s)))
    for rr in replace_rules
        sliteral = replace(sliteral, rr...)
    end
    @eval $args -> $(parse(sliteral))
end
function Sym2Function(s::Vector{SymPy.Sym}, args::Union{Symbol, Expr} = :t, replace_rules = (); fsa = false)
    vliteral = (fsa ? "Vec(" : "[") *
               join(map(x -> string(SymPy.simplify(SymPy.expand(x))), s), ", ") *
               (fsa ? ")" : "]")
    for rr in replace_rules
        vliteral = replace(vliteral, rr...)
    end
    @eval $args -> $(parse(vliteral))
end
function Sym2Function(s::Matrix{SymPy.Sym}, args::Union{Symbol, Expr} = :t, replace_rules = ())
    mliteral = "[" *
               join(mapslices(rw -> join(rw, " "), map(x -> string(SymPy.simplify(SymPy.expand(x))), s), 2), "; ") *
               "]"
    for rr in replace_rules
        mliteral = replace(mliteral, rr...)
    end
    @eval $args -> $(parse(mliteral))
end

## 2BVP Utility Functions
function LinearQuadratic2BVP{T}(A::Matrix{T}, B::Matrix{T}, c::Vector{T}, R::Matrix{T})
    n = size(A, 1)
    t = SymPy.symbols("t", real=true)
    s = SymPy.symbols("s", real=true)
    xS = [SymPy.symbols(join(["x$i" for i in 1:n], ", "), real = true)...]
    yS = [SymPy.symbols(join(["y$i" for i in 1:n], ", "), real = true)...]
    replace_rules = ((r"[xy]\d+", s -> s[1:1] * "[" * s[2:end] * "]"),
                     (r"\d+\.\d*", x -> "$T("*x*")"))

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

    LinearQuadratic2BVP(A, B, c, R,
                        Sym2Function(GinvS, :(t::$T), replace_rules),
                        Sym2Function(expAtS, :(t::$T), replace_rules),
                        Sym2Function(cdriftS, :(t::$T), replace_rules),
                        Sym2Function(costS, :(x::Vec{$n,$T}, y::Vec{$n,$T}, t::$T), replace_rules),
                        Sym2Function(dcostS, :(x::Vec{$n,$T}, y::Vec{$n,$T}, t::$T), replace_rules),
                        Sym2Function(ddcostS, :(x::Vec{$n,$T}, y::Vec{$n,$T}, t::$T), replace_rules),
                        Sym2Function(xS, :(x::Vec{$n,$T}, y::Vec{$n,$T}, t::$T, s::$T), replace_rules, fsa=true))
end

## Time-Optimal 2BVP
function topt_bisection{T}(dc, x0::Vec, x1::Vec, tm::T; tol = T(1e-3))
    # Bisection
    b = tm
    dc(x0, x1, b) < 0 && return tm
    a = tm / 100
    while dc(x0, x1, a) > 0; a /= 2; end
    m = T(0)
    cdval = T(1)
    while abs(cdval) > tol && abs(a - b) > tol
        m = (a+b)/2
        cdval = dc(x0, x1, m)
        cdval > 0 ? b = m : a = m
    end
    m
end
function topt_newton{T}(dc, ddc, x0::Vec, x1::Vec, tm::T; tol = T(1e-6))
    # Bisection / Newton's method combo
    b = tm
    dc(x0, x1, b) < 0 && return tm
    a = tm / 100
    while dc(x0, x1, a) > 0; a /= 2; end
    t = tm / 2
    cdval = dc(x0, x1, t)
    while abs(cdval) > tol && abs(a - b) > tol
        t = t - cdval / ddc(x0, x1, t)
        (t < a || t > b) && (t = (a+b)/2)
        cdval = dc(x0, x1, t)
        cdval > 0 ? b = t : a = t
    end
    t
end
function steer(L::LinearQuadratic2BVP, x0, x1, r)
    x0 == x1 && return zero(r), zero(r)
    t = isa(r, Float64) ? topt_newton(L.dcost, L.ddcost, x0, x1, r) : topt_bisection(L.dcost, x0, x1, r)
    L.cost(x0, x1, t), t
end
function steer_pairwise{T,S<:Vec}(L::LinearQuadratic2BVP{T}, V::Vector{S}, W::Vector{S}, r)
    M = length(V)
    N = length(W)
    Vmat = statevec2mat(V)
    Wmat = statevec2mat(W)
    Vbarmat = L.expAt(r)*Vmat .+ L.cdrift(r)
    Ginv = L.Ginv(r)
    BRB = L.B*inv(L.R)*L.B'

    cd = pairwise(SqMahalanobis(Ginv*BRB*Ginv), Vbarmat, Wmat)
    LHT = Ginv*(L.A*Wmat .+ L.c)

    T1 = Distances.dot_percol(Wmat, LHT)
    BLAS.gemm!('T', 'N', T(-2), Vbarmat, LHT, T(1), cd)
    broadcast!(.+, cd, cd, 2*T1')
    broadcast!(.-, cd, 1, cd)

    cands = cd .> 0
    if is(V, W)
        for i in 1:N
            cands[i,i] = false
        end
    end
    IS, JS = findn(cands)
    VTS = [steer(L, V[i], W[j], r) for (i,j) in zip(IS, JS)]
    II = find(T[v for (v,s) in VTS] .<= r)
    DS = sparse(IS[II], JS[II], T[v for (v,s) in VTS[II]], M, N)
    US = sparse(IS[II], JS[II], [DurationAndTargetControl(VTS[i][2], W[JS[i]]) for i in II], M, N)
    DS, US
end