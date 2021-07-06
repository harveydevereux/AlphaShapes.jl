using AlphaShapes, GLMakie # or using Makie if no gpu

# rotation aounrd the x axis
# anchored at [0,0,X]
function trans1(l,X)
    θ = l*π/X
    [1 0 0 0;
    0 cos(θ) -sin(θ) X*sin(θ);
    0 sin(θ) cos(θ) X-X*cos(θ);
    0 0 0 1]
end

# rotation aounrd the z axis
# anchored at [0,-Y,0]
function trans2(l,Y)
    θ = l*π/Y
    [cos(θ) -sin(θ) 0 -Y*sin(θ);
    sin(θ) cos(θ) 0 -Y+Y*cos(θ);
    0 0 1 0;
    0 0 0 1]
end

l = 1.
M = 100
N = 100

X = [x->trans1(m,l)*x for m in range(-1,stop=1,length=M)];
Y = [(trans2(n,l)*[0.0,0.0,0.0,1.0]) for n in range(-1,stop=1,length=N)];
P = []
for x in X
    for y in Y
        push!(P,x(y)[1:3])
    end
end
xs = map(x->x[1],P)
ys = map(x->x[2],P)
zs = map(x->x[3],P)

p = cat(cat(xs,ys,dims=2),zs,dims=2);

# a 3d delaunay triangulation
tets = AlphaShape(p,α=1.)
# helper function from utils
v,f = AlphaShapes.Tetrahedrons2Mesh(tets);

v = Node(v)
f = Node(f)

scene = Scene()
mesh!(scene,v,f)
scene

# play with it
GLMakie.record(scene, "makie.mp4", [1], framerate=30) do i
end
