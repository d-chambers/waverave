using WaveRave
using Test

@testset "stencil" begin
    include("test_stencil.jl")
end

@testset "control" begin
    include("test_control.jl")
end

@testset "decompose" begin
    include("test_decompose.jl")
end


@testset "utils" begin
    include("test_utils.jl")
end

@testset "2D E2@" begin
    include("test_2d_end_to_end.jl")
end