# AlphaShapes.jl
Basic implementation of alpha shapes in 2D

### Functionality:

AlphaShape is a function which returns the alpha shape triangulation from a set of 2D points

```julia
AlphaShape(X;α=nothing,search=(0.0, 10.0),MaxSteps=100)
```
Returns a list of triangles (3 by 2) arrays with vertices addressed by [n,:].

Search is a tuple passed to the optimises as the range of alpha values to search, whilst MaxSteps is the soft maximum
number of steps the optimiser takes (may take many more if each iteration is very quick)

```julia
AlphaShapeArea(A)
```

Computes the area of the alpha shape from the triangulation returned from AlphaShape
### Usage:

Given a set of 2D points, X ~ npoints by 2, compute the alpha shape (by automatically optimising alpha) using 
```Julia
A = AlphaShape(X)
```
To choose an alpha value simply pass it 
```Julia
A = AlphaShape(X,α=1.0)
```

### Uses: 
- [Triangle.jl](https://github.com/cvdlab/Triangle.jl) for Delaunay Triangulations
- [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) for alpha value optimisation

### "Normal donut" [Example](https://github.com/harveydevereux/AlphaShapes.jl/blob/master/examples/examples.jl) 


![example](https://github.com/harveydevereux/AlphaShapes.jl/blob/master/examples/Example.png)


