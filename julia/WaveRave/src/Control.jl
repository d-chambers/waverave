"""
General control functions for WaveRave
"""

using Base
include("constants.jl")


"""
Abstract type for a source. i.e. a Ricker wavelet.
"""
abstract type AbstractSource end


"""
A Ricker wavelet source.
    Source information.
    position should be a nd-array with the correct dimensions.
    time_delay controlls when the source fires in simulation time.
    frequency is in Hz.
    amplitude is the amplitude of the source.
"""
Base.@kwdef struct RickerSource <: AbstractSource
    position:: AbstractArray{Real, 1}
    time_delay:: Real = 0.0
    frequency:: Real = 10.0
    amplitude:: Real = 1.0
end


"""
A struct to hold the control parameters for the simulation

# Arguments
    p_velocity: An array of velocities in (m/s)
    density: An array of density values in (kg/m^3)
    x_values: Labels for each axis. Must have a length of the same as the
        dimension of the velocity and density arrays.
    sources: An array of sources to fire in the simulation.
    dt: time step in seconds
"""
Base.@kwdef struct WaveSimulation
    p_velocity:: OneToThreeDimRealArray
    density:: OneToThreeDimRealArray
    x_values:: Array{Array{Real, 1}}
    sources:: Array{AbstractSource, 1}
    dt:: Real
    # Optional parameters
    cfl_limit:: Real = 0.5
    nodes_per_wavelength:: Int = 10

    function WaveSimulation(velocity, density, x_values, dt)
        """
        Inner constructor for WaveControl
        """
        wc = new(velocity, density, x_values, dt)
        # validate CFL, 
        return wc
    end
 end
