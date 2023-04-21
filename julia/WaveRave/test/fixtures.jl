"""
Fixtures for testing.
"""


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
        x_values=x_values,
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
        x_values=x_values,
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
    function get_velocity_and_coords()
        x_coord = [0:dx:500]
        z_coord = [0:dx:300]
        vel = ones(Float64, length(z_coord), length(x_coord))



    end
    function get_coords()
    end
    function get_source()
    end

    coords, velocity = get_velocity_and_coords()

end