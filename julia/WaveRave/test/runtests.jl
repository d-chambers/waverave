using WaveRave
using Test

@testset "stencil" begin
    include("test_stencil.jl")
end

@testset "control" begin
    include("test_control.jl")
end
