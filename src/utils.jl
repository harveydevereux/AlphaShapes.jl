function DrawTri!(T)
        plot!([T[1,1],T[2,1]],[T[1,2],T[2,2]],label="")
        plot!([T[2,1],T[3,1]],[T[2,2],T[3,2]],label="")
        plot!([T[3,1],T[1,1]],[T[3,2],T[1,2]],label="")
end

function Ellipse(n;a=1,b=1,max=2π,θ=0,center=[0.0,0.0])
    x = [a*cos(t)*cos(θ) - b*sin(t)*sin(θ) for t in 0:(2π/n):2π] .+ center[1]
    y = [a*cos(t)*sin(θ) + b*sin(t)*cos(θ) for t in 0:(2π/n):2π] .+ center[2]
    return cat(x,y,dims=2)
end

function PackingFraction(X,ParticleArea)
    return (size(X,2)*ParticleArea) ./ Alphashapes.AlphaShapeVolume(AlphaShape(X,FindAlpha(X)))
end

"""
    Tetrahedrons2Mesh

Takes in an Array{Float64,3} of tetrahedrons, i.e an array
of (n,4,3) floats

Produces a Tuple of vertices and faces to be comprehended as
a triangle mesh, e.g Makie.mesh(vertices,faces)
"""
function Tetrahedrons2Mesh(tets::Array{Float64,3})
    faces_idxs = [[1,2,3], [1,2,4], [1,3,4], [2,4,3]]
    n = size(tets,1)
    vertices = zeros(n*4*3,3) # four faces by three vertices
    faces = zeros(n*4,3)
    for k in 1:size(tets,1)
        t = tets[k,:,:]
        for (i,f) in enumerate(faces_idxs)
            for j in 1:3
                # everyone loves flattend array indexing!
                vertices[(k-1)*4*3+(i-1)*3+j,:] = t[f[j],:]
                faces[(k-1)*4+i,j] = (k-1)*4*3+(i-1)*3+j
                # could reuse some verts to save space
            end
        end
    end
    return vertices,Int.(faces) # put these guys into: Makie.mesh(vertices,faces)
end
