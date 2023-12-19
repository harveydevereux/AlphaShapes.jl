"""
    ContiguousVoronoiCell(point_index::Int,tess,tess_inds)

For a data point (point_index) return the simplices which share this
point as a common vertex
"""
function ContiguousVoronoiCell(point_index::Int,tess,tess_inds)
    tris = []
    pos = []
    for ind in findall(x->x.==point_index,tess_inds)
        push!(tris,ind[2])
        push!(pos,ind[1])
    end
    return tess[tris,:,:],pos
end

"""
    ContiguousVoronoiCellArea(point_index::Int,tess,tess_inds)

For a data point (point_index) return the sum of the areas of the simplices
which share this point as a common vertex
"""
function ContiguousVoronoiCellArea(point_index::Int,tess,tess_inds)
    tris,inds = ContiguousVoronoiCell(point_index,tess,tess_inds)
    A = 0.0
    for i in 1:size(tris,1)
        A += AlphaShapes.SimplexVolume(tris[i,:,:])
    end
    return A
end

"""
    DTFELocalDensity(point_index::Int,tess,tess_inds)

Return the local desnity at a point in the Delaunay Tesselation Field
Estimator sense.

https://arxiv.org/pdf/0708.1441.pdf
"""
DTFELocalDensity(point_index::Int,tess,tess_inds)::Float64 = size(tess,2)/ContiguousVoronoiCellArea(point_index,tess,tess_inds)

function LawOfCosine(a::Float64,b::Float64,c::Float64)::Float64
    return acos(max(-1.0,min(1.0,(a^2 + b^2 - c^2) / (2*a*b+1e-100))))
end

function TriangleAngles(triangle::Array{Float64,2})::Array{Float64,1}
    p_a, p_b, p_c = triangle[1,:],triangle[2,:],triangle[3,:]
    a = norm(p_c .- p_b)
    b = norm(p_a .- p_c)
    c = norm(p_a .- p_b)

    A = LawOfCosine(b,c,a)
    B = LawOfCosine(a,c,b)
    C = π - (A + B)
    return [A,B,C]
end

PointWeights(triangle::Array{Float64,2})::Array{Float64,1} =  TriangleAngles(triangle) ./ 2π

"""
    WeightedDTFELocalDensity(point_index::Int,tess,tess_inds)

Return the local desnity at a point in the Weighted Delaunay Tesselation Field
Estimator sense. This additionally normalises by the 'contribution' of individual
points (portion of 2pi rad. in triangle)

 https://royalsocietypublishing.org/doi/10.1098/rsif.2021.0114
"""

function WeightedDTFELocalDensity(point_index::Int,tess,tess_inds)::Float64
    tris,pos = ContiguousVoronoiCell(point_index,tess,tess_inds)
    A = zeros(size(tris,1))
    w = zeros(size(tris,1))
    for i in 1:size(tris,1)
        w[i] = PointWeights(tris[i,:,:])[pos[i]]
        A[i] = AlphaShapes.SimplexVolume(tris[i,:,:])
    end
    return sum(w./A)/2.0*π
end
