using AlphaShapes
using Test

@testset "Test Volume" begin
   p = [0 1;1 0;1 1]
   @test AlphaShapes.SimplexVolume(p) == 0.5
end
