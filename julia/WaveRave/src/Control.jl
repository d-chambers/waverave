"""
General control functions for WaveRave
"""

using Base
using Statistics

include("Sources.jl")


"""
A struct to hold the control parameters for the simulation

# Arguments
    dx: dx values along each dimension in (m)
    dt: time step in seconds
    p_velocity: An array of velocities in (m/s)
    density: An array of density values in (kg/m^3)
    sources: An array of sources to fire in the simulation.
    cfl_limit: The CFL limit to enforce. 0.5 is recommended.
    nodes_per_wavelength: The number of nodes per wavelength to enforce.
        The slowest velocity in the simulation and highest frequency source
        will be used to calculate the wavelength.
"""
Base.@kwdef struct WaveSimulation
    dt :: Real
    x_values :: AbstractArray{AbstractArray}
    p_velocity :: AbstractArray
    density :: AbstractArray
    sources :: Array{AbstractSource, 1}
    # Optional parameters
    cfl_limit:: Real = 0.5
    nodes_per_wavelength:: Int = 10
    processors:: Int = 1
    division:: String = "grid"
 end


# --- Validation functions ---

"""
    A function to run a quite of validation on the simulation.
        Currently these include:
        1. Ensure all arrays have consistent dimensions.
        2. Ensure spatial values have the same shape as the property arrays.
        3. Ensure that the spatial values are monotonic and evenly sampled.
        4. Check CFL limit is not exceeded. 
        5. Check that nodes_per_wavelegnth is met.
"""
 function validate_simulation(simulation:: WaveSimulation)
    check_array_dimensions(simulation)
    check_xvalue_lengths(simulation)
    check_xvalues(simulation)
    check_cfl(simulation)
    check_wavelength_resolution(simulation)
    return true
 end


"""
    Check that all arrays have the same dimensions.
"""
function check_array_dimensions(simulation:: WaveSimulation)
    if !all(size(simulation.p_velocity) .== size(simulation.density))
        error("p_velocity and density arrays must have the same dimensions")
    end
end


"""
Check that the spatial values have consistent shapes with the property arrays.
"""
function check_xvalue_lengths(simulation:: WaveSimulation)
    grid_shape = size(simulation.p_velocity)
    for (ind, x) in enumerate(simulation.x_values)
        if size(x)[1] != grid_shape[ind]
            error("x_value lengths not consistent with grid shape in dimension $ind")
        end
    end
end


"""
Check that the spatial values are monotonic and evenly sampled.
"""
function check_xvalues(simulation:: WaveSimulation)
    for (ind, x) in enumerate(simulation.x_values)
        diffs = diff(x)
        if !all(diffs .> 0)
            error("x_values must be monotonic!, they are not in dimension $ind")
        end
        if !all(diffs .≈ median(diffs))
            error("x_values must be evenly sampled!, they are not in dimension $ind")
        end
    end
end
   


 """
    Check the CFL limit is not exceeded.
 """
function check_cfl(simulation:: WaveSimulation)
    dx_min = minimum([minimum(diff(x)) for x in simulation.x_values])
    max_velocity = maximum(simulation.p_velocity)
    dt = simulation.dt
    if dx_min/dt < max_velocity * simulation.cfl_limit
        info = "$(dx_min/dt) < $max_velocity * $(simulation.cfl_limit)"
        error("CFL limit exceeded! $info")
    end
    
end


"""
Check wavelength resolution is met.
"""
function check_wavelength_resolution(simulation:: WaveSimulation)
    # Get the slowest velocity and highest frequency source
    slowest_velocity = minimum(simulation.p_velocity)
    highest_frequency = maximum([source.frequency for source in simulation.sources])
    wavelength = slowest_velocity / highest_frequency
    dx_min = minimum([minimum(diff(x)) for x in simulation.x_values])
    if wavelength / dx_min < simulation.nodes_per_wavelength
        info = "$wavelength / $dx_min < $(simulation.nodes_per_wavelength)"
        error("Wavelength resolution not met! $info")
    end
end
