using MotionPlanning

const ISRR_2H = Compound2D(
    Shape2D{Float64}[
        Box2D([0., 0.16], [0.36, 0.5]),
        Box2D([0.4, 0.5], [0.19, 0.35]),
        Box2D([0.22, 0.46], [0.57, 0.75]),
        Box2D([0.75, 1.], [0.64, 0.77]),
        Box2D([0.22, 0.8], [0.34, 0.51])
    ]
)

const TRI_BALLS = Compound2D(
    Shape2D{Float64}[
        Polygon([(.3,.3), (.7,.3), (.5,.65)]),
        Circle([.3,.3], .15),
        Circle([.7,.3], .15),
        Circle([.5,.65], .15)
    ]
)

const ISRR_POLY = Compound2D(
    Shape2D{Float64}[
        Polygon([(0., .25), (.27, .28), (.17, .4), (0., .4)]),
        Polygon([(.5, .2), (.2, .5), (.25, .7), (.4, .8), (.6, .8), (.7, .5)]),
        Polygon([(.55, .2), (.75, .5), (.85, .5), (.85, .2)]),
        # Polygon([(.3, .6), (.15, .85), (.4, .6)]),
        Circle([.9, .65], .1)
    ]
)

const ISRR_POLY_WITH_SPIKE = Compound2D(
    Shape2D{Float64}[
        Polygon([(0., .25), (.27, .28), (.17, .4), (0., .4)]),
        Polygon([(.5, .2), (.2, .5), (.25, .7), (.4, .8), (.6, .8), (.7, .5)]),
        Polygon([(.55, .2), (.75, .5), (.85, .5), (.85, .2)]),
        Polygon([(.3, .6), (.15, .85), (.4, .6)]),
        Circle([.9, .65], .1)
    ]
)

const EMPTY_2D = Compound2D(Shape2D{Float64}[])

nothing