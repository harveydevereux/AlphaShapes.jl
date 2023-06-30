using AlphaShapes
using Test
using LinearAlgebra
include("../src/utils.jl")

@testset "Test Alpha Shape 2D" begin
   N = 1000
   r1 = 0.75
   r2 = 2.0
   X = randn(2, N)
   ind = [((norm(X[:,i])>r1) + (norm(X[:,i])<r2)) == 2 for i in axes(X,2)]
   X = X[:,ind]
   AlphaShape(X;α=nothing,search=(0.0, 10.0),MaxSteps=100)
   t = (@timed AlphaShape(X;α=nothing,search=(0.0, 10.0),MaxSteps=100))[2]
   @info "Runtime for 1000 2D points and 100 optimisation steps", t
end
