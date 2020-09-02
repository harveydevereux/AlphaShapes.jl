"""
    ContiguousVoronoiCell(point_index::Int,tess,tess_inds)

For a data point (point_index) return the simplices which share this
point as a common vertex
"""
function ContiguousVoronoiCell(point_index::Int,tess,tess_inds)
    tris = [ind[2] for ind in findall(x->x.==point_index,tess_inds)]
    return tess[tris,:,:]
end

"""
    ContiguousVoronoiCellArea(point_index::Int,tess,tess_inds)

For a data point (point_index) return the sum of the areas of the simplices
which share this point as a common vertex
"""
function ContiguousVoronoiCellArea(point_index::Int,tess,tess_inds)
    tris = ContiguousVoronoiCell(point_index,tess,tess_inds)
    A = 0.0
    for i in 1:size(tris,1)
        A += AlphaShapes.SimplexVolume(tris[i,:,:])
    end
    return A
end

"""
    LocalDensity(point_index::Int,tess,tess_inds)

Return the local desnity at a point in the Delaunay Tesselation Field
Estimator sense.

https://arxiv.org/pdf/0708.1441.pdf
"""
LocalDensity(point_index::Int,tess,tess_inds) = size(tess,2)/CellArea(point_index,tess,tess_inds)
