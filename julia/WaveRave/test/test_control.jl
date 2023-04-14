"""
Tests for the control module.
"""

using WaveRave
using Test
using Debugger

"""
A function to get a simple 1D control structure for testing.
"""
function get_simple_1d_control()
    # Create a simple 1D control structure
    p_velocity = fill(5_000., 1_000)
    density = fill(1000., sizeof(p_velocity))
    x_values = [[0:0.1:100]]
    dt = 0.1
    return WaveSimulation(p_velocity, density, x_values, dt)
end
     


@testset "test 1d" begin
    @testset "simple init" begin
        @enter wc = get_simple_1d_control()
        @test isa(wc, WaveSimulation)
    end
end
  