"""
Module for testing utils.
"""

using Test
using Debugger

using WaveRave: find_division_shape


@testset "test 1D" begin
    # happy paths
    @test find_division_shape(10, (100,)) == ([10], [10])
    @test find_division_shape(1, (100,)) == ([1], [100])
    # ensure error thrown
    @test_throws "cannot be divided" find_division_shape(11, (100,))
    # ensure 1 works

end

@testset "test 2D" begin
    # happy paths
    @enter find_division_shape(4, (100,100))
    @enter find_division_shape(4, (100,100)) == ([2, 2], [50, 50])
    # @test find_division_shape(2, (100,100)) == [100, 50]
    # @test find_division_shape(6, (300,200)) == [100, 100]
    # # ensure error thrown
    # @test_throws "cannot be divided" find_division_shape(11, (100,))
end
