"""
Fixtures for testing.
"""


"""
A function to get a simple 1D control structure for testing.
"""
function get_simple_1d_control()
    # Create a simple 1D control structure
    p_velocity = fill(5_000., 1_000)
    density = fill(1000., size(p_velocity))
    dt = 0.001
    x_values = [Vector(0:p_velocity[1]*dt*2.:10*(length(p_velocity)-1))]
    source = ExplosiveRickerSource(location=[50.0])
    return WaveSimulation(
        dt=dt,
        x_values=x_values,
        p_velocity=p_velocity, 
        density=density, 
        sources=[source], 
    )
end


"""
A function to get a controll structure which doesn't meet cfl limit.
"""
function get_simple_1d_bad_clf()
    # Create a simple 1D control structure
    p_velocity = fill(5_000., 1_000)
    density = fill(1000., size(p_velocity))
    dt = 0.1
    x_values = [Vector(0:1:1*(length(p_velocity)-1))]
    source = ExplosiveRickerSource(location=[50.0])
    return WaveSimulation(
        dt=dt,
        x_values=x_values,
        p_velocity=p_velocity, 
        density=density, 
        sources=[source], 
    )
end
