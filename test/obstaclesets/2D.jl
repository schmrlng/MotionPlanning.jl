using MotionPlanning

ISRR_2H{T<:AbstractFloat}(::Type{T} = Float64) = Compound2D(
    Shape2D{T}[
        Box2D([T(.0), T(.16)], [T(.36), T(.5)]),
        Box2D([T(.4), T(.5)], [T(.19), T(.35)]),
        Box2D([T(.22), T(.46)], [T(.57), T(.75)]),
        Box2D([T(.75), T(1.)], [T(.64), T(.77)]),
        Box2D([T(.22), T(.8)], [T(.34), T(.51)])
    ]
)

TRI_BALLS{T<:AbstractFloat}(::Type{T} = Float64) = Compound2D(
    Shape2D{T}[
        Polygon([(T(.3), T(.3)), (T(.7), T(.3)), (T(.5), T(.65))]),
        Circle([T(.3), T(.3)], T(.15)),
        Circle([T(.7), T(.3)], T(.15)),
        Circle([T(.5), T(.65)], T(.15))
    ]
)

ISRR_POLY{T<:AbstractFloat}(::Type{T} = Float64) = Compound2D(
    Shape2D{T}[
        Polygon([(T(.0), T(.25)), (T(.27), T(.28)), (T(.17), T(.4)), (T(.0), T(.4))]),
        Polygon([(T(.5), T(.2)), (T(.2), T(.5)), (T(.25), T(.7)), (T(.4), T(.8)), (T(.6), T(.8)), (T(.7), T(.5))]),
        Polygon([(T(.55), T(.2)), (T(.75), T(.5)), (T(.85), T(.5)), (T(.85), T(.2))]),
        # Polygon([(T(.3), T(.6)), (T(.15), T(.85)), (T(.4), T(.6))]),
        Circle([T(.9), T(.65)], T(.1))
    ]
)

ISRR_POLY_WITH_SPIKE{T<:AbstractFloat}(::Type{T} = Float64) = Compound2D(
    Shape2D{T}[
        Polygon([(T(.0), T(.25)), (T(.27), T(.28)), (T(.17), T(.4)), (T(.0), T(.4))]),
        Polygon([(T(.5), T(.2)), (T(.2), T(.5)), (T(.25), T(.7)), (T(.4), T(.8)), (T(.6), T(.8)), (T(.7), T(.5))]),
        Polygon([(T(.55), T(.2)), (T(.75), T(.5)), (T(.85), T(.5)), (T(.85), T(.2))]),
        Polygon([(T(.3), T(.6)), (T(.15), T(.85)), (T(.4), T(.6))]),
        Circle([T(.9), T(.65)], T(.1))
    ]
)

EMPTY_2D{T<:AbstractFloat}(::Type{T} = Float64) = Compound2D(Shape2D{T}[])

nothing