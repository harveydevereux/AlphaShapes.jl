# AlphaShapes.jl
Basic implementation of alpha shapes in 2D
| Travis | AppVeyor | CodeCov |
|-------|----------|---------|
| [![Build Status](https://travis-ci.com/harveydevereux/AlphaShapes.jl.svg?branch=master)](https://travis-ci.com/harveydevereux/AlphaShapes.jl)|[![Build status](https://ci.appveyor.com/api/projects/status/r3x53cnn3n7352tf?svg=true)](https://ci.appveyor.com/project/harveydevereux/AlphaShapes-jl) |[![codecov](https://codecov.io/gh/harveydevereux/DNAText.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/harveydevereux/AlphaShapes.jl) |

### Functionality:

AlphaShape is a function which returns the alpha shape from a set of N dimensional points

```julia
AlphaShape(X;α=nothing,search=(0.0, 10.0),MaxSteps=100)
```
Returns a list of triangles (3 by N+1) arrays with vertices addressed by [v,:].

Search is a tuple passed to the optimiser as the range of alpha values to search, whilst MaxSteps is the soft maximum
number of steps the optimiser takes (may take many more if each iteration is very quick)

```julia
AlphaShapeArea(A)
```

Computes the area of the alpha shape from the triangulation returned from AlphaShape
### Usage:

Given a set of N-D points, X ~ npoints by N, compute the alpha shape (by automatically optimising alpha) using 
```Julia
A = AlphaShape(X)
```
To choose an alpha value simply pass it 
```Julia
A = AlphaShape(X,α=1.0)
```

### Uses: 
- [MiniQhull.jl](https://github.com/gridap/MiniQhull.jl) for Delaunay Triangulations in any dimension
- [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) for alpha value optimisation

### "Normal donut" [Example](https://github.com/harveydevereux/AlphaShapes.jl/blob/master/examples/examples.jl) 


![example](https://github.com/harveydevereux/AlphaShapes.jl/blob/master/examples/Example.png)


