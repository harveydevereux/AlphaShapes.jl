# AlphaShapes.jl Documentation

## Overview of AlphaShapes.jl

### Alpha Shapes

An $\alpha$-shape is closely related to the Delaunay triangulation and Convex Hull
of a set of points (in any dimension).

The $\alpha$-shape is constructed for a parameter $\alpha$ by testing each
Delaunay simplex to determine if it's circumsphere has a radius less than
$\alpha$. If it does then the simplex is kept, otherwise the simplex is
removed. This modified triangulation is the $\alpha$-shape.

Choosing $\alpha$ is done in many ways. AlphaShapes.jl optimises $\alpha$
to find the triangulation with the smallest volume which includes all input
points as a vertex to at least on simplex.

[A mathematical reference](https://graphics.stanford.edu/courses/cs268-11-spring/handouts/AlphaShapes/as_fisher.pdf)

### Functionality

AlphaShapes.jl represents an $\alpha$-shape as a function of an array,
```Array{Float64,2}```, of points. So all user functionality comes through this.

```@docs
AlphaShape
```

### Utils

The volume of one simplex can be found using

```@docs
AlphaShapes.SimplexVolume
```

For the whole $\alpha$-shape the volume can be found using

```@docs
AlphaShapes.AlphaShapeVolume
```

The Delaunay triangulation is used to build the alpha shapes so there
is a helper method to build one for any set of points (MiniQhull.jl)

```@docs
AlphaShapes.GetDelaunayTriangulation
```

### Under the hood

Arbitrary dimensions is handled by using the Cayley-Menger matrix for simplex
volume and simplex circumsphere calculation

```@docs
AlphaShapes.CayleyMenger
AlphaShapes.SimplexCircumSphere
```


## contents

```@contents
```

## Index

```@index
```
