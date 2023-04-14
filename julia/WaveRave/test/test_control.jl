"""
Tests for the control module.
"""

using WaveRave
using Test
using Debugger

# -------- Fixtures --------

"""
A function to get a simple 1D control structure for testing.
"""
function get_simple_1d_control()
    # Create a simple 1D control structure
    p_velocity = fill(5_000., 1_000)
    density = fill(1000., size(p_velocity))
    dt = 0.001
    x_values = [Vector(0:p_velocity[1]*dt*2.:10*(length(p_velocity)-1))]
    source = ExplosiveRickerSource(position=[50.0])
    return WaveSimulation(
        dt=dt,
        x_values=x_values,
        p_velocity=p_velocity, 
        density=density, 
        sources=[source], 
    )
end


"""
A function to get a controll structure which doesn't meet clfl limit.
"""
function get_simple_1d_bad_clf()
    # Create a simple 1D control structure
    p_velocity = fill(5_000., 1_000)
    density = fill(1000., size(p_velocity))
    dt = 0.1
    x_values = [Vector(0:1:1*(length(p_velocity)-1))]
    source = ExplosiveRickerSource(position=[50.0])
    return WaveSimulation(
        dt=dt,
        x_values=x_values,
        p_velocity=p_velocity, 
        density=density, 
        sources=[source], 
    )
end


# ---- Tests -----


@testset "test 1d" begin
    @testset "simple init" begin
        wc = get_simple_1d_control()
        @test isa(wc, WaveSimulation)
    end
    @testset "test CFL limit" begin
        # case that should pass
        wc = get_simple_1d_control()
        @test validate_simulation(wc)
        # case that should fail
        wc = get_simple_1d_bad_clf()
        @test_throws "CFL limit exceeded!" validate_simulation(wc)
    end
end
  
