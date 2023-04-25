"""
Module for testing utils.
"""

using Test
using Debugger

import WaveRave
using WaveRave: find_division_shape


@testset "test 1D find_division_shape" begin
    # happy paths
    @test find_division_shape(10, (100,)) == ([10], [10])
    @test find_division_shape(1, (100,)) == ([1], [100])
    # ensure error thrown
    @test_throws "cannot be divided" find_division_shape(11, (100,))
    # ensure 1 works

end


@testset "test 2D find_division_shape" begin
    # happy paths
    @test find_division_shape(4, (100,100)) == ([2, 2], [50, 50])
    @test find_division_shape(2, (100,100)) == ([2, 1], [50, 100])
    @test find_division_shape(6, (300,200)) == ([3, 2], [100, 100])
    # ensure error thrown
    @test_throws "cannot be divided" find_division_shape(7, (100, 100))
end


@testset "test make homogeneous model" begin
    mod = WaveRave.make_homogeneous_model((300, 500), 5, 1_000)
    @test size(mod) == (60, 100)
end
