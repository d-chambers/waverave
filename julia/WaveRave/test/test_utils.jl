"""
Module for testing utils.
"""

using Test
using Debugger

using WaveRave: find_division_shape

@testset "test divisor shape" begin
    @testset "test 1D" begin
        out = find_division_shape(10, (100,))
        @test out == [10]
    end
end
