using AlphaShapes
using Test

@testset "Test Local Density" begin
    # four corners
    X = [0.0 0.0 1.0 1.0; 0.0 1.0 1.0 0.0]
    tess, inds = AlphaShapes.GetDelaunayTriangulation(X,true)
    density = [AlphaShapes.WeightedDTFELocalDensity(i,tess,inds) for i in 1:size(X,2)]

    @test density == [
        0.39269908169872425,
        0.39269908169872403,
        0.39269908169872425,
        0.39269908169872403
    ]

    density = [AlphaShapes.DTFELocalDensity(i,tess,inds) for i in 1:size(X,2)]

    @test density == [
        3.0,
        1.5,
        3.0,
        1.5       
    ]
end