const Shape2D = SeparatingAxisTheorem2D.Shape2D

intersecting(x::Shape2D, y         , args...) = SeparatingAxisTheorem2D.intersecting(x, y, args...)
# intersecting(x         , y::Shape2D, args...) = SeparatingAxisTheorem2D.intersecting(x, y, args...)
# intersecting(x::Shape2D, y::Shape2D, args...) = SeparatingAxisTheorem2D.intersecting(x, y, args...)

sweep_intersecting(x::Shape2D, y         , args...) = SeparatingAxisTheorem2D.sweep_intersecting(x, y, args...)
# sweep_intersecting(x         , y::Shape2D, args...) = SeparatingAxisTheorem2D.sweep_intersecting(x, y, args...)
# sweep_intersecting(x::Shape2D, y::Shape2D, args...) = SeparatingAxisTheorem2D.sweep_intersecting(x, y, args...)
