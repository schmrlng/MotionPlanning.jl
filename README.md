# MotionPlanning.jl

[![Build Status](https://travis-ci.org/schmrlng/MotionPlanning.jl.svg?branch=master)](https://travis-ci.org/schmrlng/MotionPlanning.jl)

A Julia package containing motion planning code related to the following papers on robotic motion planning:
- [Fast Marching Tree: a Fast Marching Sampling-Based Method for Optimal Motion Planning in Many Dimensions](http://arxiv.org/abs/1306.3532)
- [Optimal Sampling-Based Motion Planning under Differential Constraints: the Driftless Case](http://arxiv.org/abs/1403.2483)
- [Optimal Sampling-Based Motion Planning under Differential Constraints: the Drift Case with Linear Affine Dynamics](http://arxiv.org/abs/1405.7421)
- [Monte Carlo Motion Planning for Robot Trajectory Optimization Under Uncertainty](http://arxiv.org/abs/1504.08053)
- [Evaluating Trajectory Collision Probability through Adaptive Importance Sampling for Safe Motion Planning](https://arxiv.org/abs/1609.05399)

Install with
```julia
Pkg.clone("https://github.com/schmrlng/MotionPlanning.jl.git")
```
at the ```julia>``` prompt (and let me know if it actually works!).

## Documentation
Doesn't exist. At least for now. Check out some [basic usage examples](http://nbviewer.ipython.org/github/schmrlng/MotionPlanning.jl/blob/master/doc/MotionPlanning.ipynb) and peruse the source, or message me if you have any particular questions about the code.

## Disclaimers
- This code is so alpha that if it was particle, it never would have made it past the paper planning phase.
- This code is so alpha that if it was a gorilla, even King Kong would be intimidated by its immense load time.
- This code is so alfa that if it was a car, Jeremy Clarkson would revel in its shoddy engineering.
- This code is so alfalfa... you get the picture.
In all seriousness, lots of things in this package are likely to change in the immediate future as I port in old code, reorganize things, and address some of the 50+ "TODO"s scattered throughout the code.
