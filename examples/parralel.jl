using Distributed

using AlphaShapes

T = randn(1000,3)

@time AlphaShape(T)

rmprocs(workers())
addprocs(6)
@info workers()
@everywhere using BlackBoxOptim
@everywhere using AlphaShapes

# should be faster!
@time AlphaShape(T)

rmprocs(workers())
