using Test
using WaveRave


@testset "test centered difference first derivative" begin
    # test for 3 point stencil
    w = get_stencil(-1:1, 0, 1)
    @test w ≈ [-1/2, 0, 1/2]
    # test for 5 point stencil 
    w = get_stencil(-2:2, 0, 1)
    @test w ≈ ([1, -8, 0, 8, -1] ./ 12)
end


@testset "test centered difference second derivative" begin
    w = get_stencil(-1:1, 0, 2)
    @test w ≈ [1, -2, 1]
end


@testset "test 0th deriviative" begin
    w = get_stencil(-1:1, 0, 0)
    @test w ≈ [0, 1, 0]
    w = get_stencil(0:3, 0, 0)
    @test w ≈ [1, 0, 0, 0]
end