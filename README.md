# AlphaShapes.jl Dev

Goals: (1) 3D quality of life, visualisations and animations with Makie, (2) Parralelism (alpha value determination)

Progress notes:

- 3D usecases, [e.g](https://github.com/harveydevereux/AlphaShapes.jl/blob/dev/examples/torus-makie.jl) 
  - seems to be an increase in errors in CM inversion
  - alpha determination seems more difficult


![torus](https://github.com/harveydevereux/AlphaShapes.jl/blob/dev/examples/torus.png)


- Parallelism
  - Not at all seamlessly integrated, but can be used if the user is willing to follow [this example](https://github.com/harveydevereux/AlphaShapes.jl/blob/dev/examples/parralel.jl) 
