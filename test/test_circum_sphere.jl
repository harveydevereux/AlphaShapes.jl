using AlphaShapes
using Test

tol = 1e-6

@testset "Test Circum Sphere" begin
   p = [0. 1.;1. 0.;1. 1.]
   c,r = AlphaShapes.SimplexCircumSphere(p)
   @test (r .- 0.7071067811865476) < tol
   @test (c[1] .- 0.5) < tol
   @test (c[1].-0.5) < tol
end
