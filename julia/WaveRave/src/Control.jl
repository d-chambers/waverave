"""
General control functions for WaveRave
"""

using Base
using Statistics

include("Sources.jl")
include("Receivers.jl")


"""
    Get dt based on coordinates.

This is a little quick and dirty; we should probably match
dx with specific velocities. 
"""
function get_dt(coords, velocity, cfl_limit)
    max_vel = maximum(velocity)
    min_dx = minimum([minimum(x[2:end] - x[1:end - 1]) for x in coords])
    return (min_dx / max_vel) * cfl_limit
end


"""
    A struct to hold the control parameters for the simulation

# Arguments
    coords: coordinates along each dimension. Should be evenly sampled.
    p_velocity: An array of velocities in (m/s)
    sources: An array of sources to fire in the simulation.
    dt: time step in seconds, if None calculate based on CFL.
    receivers: An array of receivers for recording wavefields.
    cfl_limit: The CFL limit to enforce. 0.5 is recommended.
    nodes_per_wavelength: The number of nodes per wavelength to enforce.
        The slowest velocity in the simulation and highest frequency source
        will be used to calculate the wavelength.
    distribution_strategy: The strategy to use for distributing the simulation
        across the processors. options are: 
        :grid - The domian is decomposed into evenly sized block subdomains.
        :stip - The domain is decomposed into evenly sized strip subdomains.
    distributed_shape: The shape of each distributed chunk.
    space_order: The number of points right of the centeral point for the 
        derivative estimation in space. Centeral is always used.
    time_order: The number of points right of the centeral point for the 
        derivative estimation in space. Centeral is always used.
"""
Base.@kwdef mutable struct WaveSimulation
    coords::AbstractArray{AbstractArray}
    p_velocity::AbstractArray
    sources::Array{AbstractSource,1}
    # Optional parameters
    receivers::Array{Receiver, 1} = []
    cfl_limit::Real = 0.5
    nodes_per_wavelength::Int = 10
    distribution_strategy:: Symbol = :grid
    distributed_shape:: Union{Tuple, Nothing} = nothing
    space_order::Int = 5
    time_order::Int = 1
    dt::Real = get_dt(coords, p_velocity, cfl_limit)
end


function (wavesim::WaveSimulation)()
    validate_simulation(wavesim)
    return wavesim
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
function validate_simulation(simulation::WaveSimulation)
    # check_array_dimensions(simulation)
    check_xvalue_lengths(simulation)
    check_xvalues(simulation)
    check_cfl(simulation)
    check_wavelength_resolution(simulation)
    return true
end


# """
#     Check that all arrays have the same dimensions.
# """
# function check_array_dimensions(simulation::WaveSimulation)
#     if !all(size(simulation.p_velocity) .== size(simulation.density))
#         error("p_velocity and density arrays must have the same dimensions")
#     end
# end


"""
Check that the spatial values have consistent shapes with the property arrays.
"""
function check_xvalue_lengths(simulation::WaveSimulation)
    grid_shape = size(simulation.p_velocity)
    for (ind, x) in enumerate(simulation.coords)
        if size(x)[1] != grid_shape[ind]
            error("x_value lengths not consistent with grid shape in dimension $ind")
        end
    end
end


"""
Check that the spatial values are monotonic and evenly sampled.
"""
function check_xvalues(simulation::WaveSimulation)
    for (ind, x) in enumerate(simulation.coords)
        diffs = diff(x)
        if !all(diffs .> 0)
            error("coords must be monotonic!, they are not in dimension $ind")
        end
        if !all(diffs .≈ median(diffs))
            error("coords must be evenly sampled!, they are not in dimension $ind")
        end
    end
end



"""
   Check the CFL limit is not exceeded.
"""
function check_cfl(simulation::WaveSimulation)
    dx_min = minimum([minimum(diff(x)) for x in simulation.coords])
    max_velocity = maximum(simulation.p_velocity)
    dt = simulation.dt
    if dx_min / dt < max_velocity * simulation.cfl_limit
        info = "$(dx_min/dt) < $max_velocity * $(simulation.cfl_limit)"
        error("CFL limit exceeded! $info")
    end

end


"""
    Check wavelength resolution is met.
"""
function check_wavelength_resolution(simulation::WaveSimulation)
    # Get the slowest velocity and highest frequency source
    slowest_velocity = minimum(simulation.p_velocity)
    highest_frequency = maximum([source.frequency for source in simulation.sources])
    wavelength = slowest_velocity / highest_frequency
    dx_min = minimum([minimum(diff(x)) for x in simulation.coords])
    if wavelength / dx_min < simulation.nodes_per_wavelength
        info = "$wavelength / $dx_min < $(simulation.nodes_per_wavelength)"
        error("Wavelength resolution not met! $info")
    end
end


"""
    Local control structure (each rank gets one!)

See Docs on Simulation for descriptions.
"""
Base.@kwdef struct LocalGrid
    coords::AbstractArray{AbstractArray}
    p_velocity::AbstractArray
    origin_inds::AbstractArray{StepRange}=[]
    sources::Array{AbstractSource, 1}
    receivers::Array{Receiver, 1} = []
    global_index::AbstractArray = []
end

