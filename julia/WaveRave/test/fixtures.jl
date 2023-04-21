"""
Fixtures for testing.
"""

using WaveRave


"""
A function to get a simple 1D control structure for testing.
"""
function get_simple_1d_control()
    # Create a simple 1D control structure
    p_velocity = fill(5_000., 1_000)
    dt = 0.001
    x_values = [Vector(0:p_velocity[1]*dt*2.:10*(length(p_velocity)-1))]
    source = ExplosiveRickerSource(location=[50.0])
    return WaveSimulation(
        dt=dt,
        coords=x_values,
        p_velocity=p_velocity,  
        sources=[source], 
    )
end


"""
A function to get a controll structure which doesn't meet cfl limit.
"""
function get_simple_1d_bad_clf()
    # Create a simple 1D control structure
    p_velocity = fill(5_000., 1_000)
    dt = 0.1
    x_values = [Vector(0:1:1*(length(p_velocity)-1))]
    source = ExplosiveRickerSource(location=[50.0])
    return WaveSimulation(
        dt=dt,
        coords=x_values,
        p_velocity=p_velocity, 
        sources=[source], 
    )
end


"""
    Create the end2end control for simulation case.

This is a very specific case used for the Math 540 term
project.
"""
function get_2d_test_control(dx=1)
    function get_velocity_and_coords(fast_vel=5000.0, slow_vel=3500.0)
        x_coord = collect([0:dx:499]...) .+ 0.5
        z_coord = collect([0:dx:299]...) .+ 0.5

        in_slow_region = @. (z_coord > 135) & (z_coord < 165)
        vel = fill(fast_vel, length(z_coord), length(x_coord))
        vel[in_slow_region, :] .= slow_vel
        # set layers

        return [z_coord, x_coord], vel



    end

    function get_source()
        source = ExplosiveRickerSource(
            location=[150, 250],
            frequency=20.,
            time_delay=1.0,

        )
        return [source]
    end

    coords, velocity = get_velocity_and_coords()
    sources = get_source()

    out = WaveSimulation(
        coords=coords,
        p_velocity=velocity,
        sources=sources,
    )
    return out


end