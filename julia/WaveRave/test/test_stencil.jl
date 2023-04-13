using Test
using WaveRave


@testset "test centered difference first derivative" begin
    w = get_stencil(-1:1, 0, 1)
    @test w ≈ [-1/2, 0, 1/2]
end


@testset "test centered difference second derivative" begin
    w = get_stencil(-1:1, 0, 2)
    @test w ≈ [-1/2, 0, 1/2]
end
