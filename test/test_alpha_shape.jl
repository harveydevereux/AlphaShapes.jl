using AlphaShapes
using Test
using LinearAlgebra
include("../src/utils.jl")

@testset "Test Alpha Shape 2D" begin
   N = 1000
   r1 = 0.75
   r2 = 2.0
   X = randn(N,2)
   ind = [((norm(X[i,:])>r1) + (norm(X[i,:])<r2)) == 2 for i in 1:size(X,1)]
   X = X[ind.==true,:]
   AlphaShape(X;α=nothing,search=(0.0, 10.0),MaxSteps=100)
   t = (@timed AlphaShape(X;α=nothing,search=(0.0, 10.0),MaxSteps=100))[2]
   @info "Runtime for 1000 2D points and 100 optimisation steps", t
end
