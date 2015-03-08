using MotionPlanning

ISRR_2H = Compound2D(
    Shape2D{Float64}[
        Box2D([0., 0.16], [0.36, 0.5]),
        Box2D([0.4, 0.5], [0.19, 0.35]),
        Box2D([0.22, 0.46], [0.57, 0.75]),
        Box2D([0.75, 1.], [0.64, 0.77]),
        Box2D([0.22, 0.8], [0.34, 0.51])
    ]
)

TRI_BALLS = Compound2D(
    Shape2D{Float64}[
        Polygon([(.3,.3), (.7,.3), (.5,.65)]),
        Circle([.3,.3], .15),
        Circle([.7,.3], .15),
        Circle([.5,.65], .15)
    ]
)

nothing