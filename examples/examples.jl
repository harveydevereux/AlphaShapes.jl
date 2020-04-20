include("../src/AlphaShapes.jl")
include("../src/utils.jl")
using Plots
using Main.AlphaShapes
using LinearAlgebra

N = 1000
r1 = 0.75
r2 = 2.0
X = randn(N,2)
ind = [((norm(X[i,:])>r1) + (norm(X[i,:])<r2)) == 2 for i in 1:size(X,1)]
X = X[ind.==true,:]
# normally distributed points with a whole in the middle
p1 = scatter(X[:,1],X[:,2],label="",aspect_ratio=:equal)

T = basic_triangulation_vertices(X);
p1 = scatter(X[:,1],X[:,2],label="",aspect_ratio=:equal)
[DrawTri!(T[i]) for i in 1:size(T,1)]
plot!(title="Delaunay Triangulation")

A = AlphaShape(X);
p2 = scatter(X[:,1],X[:,2],label="",aspect_ratio=:equal)
[DrawTri!(A[i]) for i in 1:size(A,1)]
plot!(title="Optimum Alpha Shape Triangulation")

p = plot(p1,p2,size=(900,300))
savefig(p,"Example.pdf")
