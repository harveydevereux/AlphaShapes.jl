using AlphaShapes, Plots, LinearAlgebra, Random

function DrawTri!(T)
    plot!([T[1,1],T[1,2]],[T[2,1],T[2,2]],label="")
    plot!([T[1,2],T[1,3]],[T[2,2],T[2,3]],label="")
    plot!([T[1,3],T[1,1]],[T[2,3],T[2,1]],label="")
end

Random.seed!(31415926535897)

N = 10000
r1 = 0.75
r2 = 2.0

X = randn(2,N)
ind = [((norm(X[:,i])>r1) + (norm(X[:,i])<r2)) == 2 for i in axes(X,2)]
X = X[:,ind.==true]
# normally distributed points with a whole in the middle
p1 = scatter(X[1,:],X[2,:],label="",aspect_ratio=:equal)

T = AlphaShapes.GetDelaunayTriangulation(X);
p1 = scatter(X[1,:],X[2,:],label="",aspect_ratio=:equal)
[DrawTri!(T[:,:,i]) for i in axes(T,3)]
plot!(title="Delaunay Triangulation")

A = AlphaShape(X);
p2 = scatter(X[1,:],X[2,:],label="",aspect_ratio=:equal)
[DrawTri!(A[:,:,i]) for i in axes(A,3)]
plot!(title="Optimum Alpha Shape Triangulation")

p = plot(p1,p2,size=(900,300))
savefig(p,"Example.png")
