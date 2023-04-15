"""
Code for solving the acoustic wave equation.
"""

include("Control.jl")
include("Sources.jl")


"""
    A worker for running a subdomain of the simulation.
"""
Base.@kwdef mutable struct Worker
    id:: Int
    density:: AbstractArray
    p_velocity:: AbstractArray
    x_values:: AbstractArray
    sources:: Array{AbstractSource,1}
    pressure_0:: AbstractArray = zeros(size(density))
    pressure_1:: AbstractArray = zeros(size(density)) 
end



"""
    Generate a worker for a given ID based on the 
"""
function Worker(simulation::WaveSimulation, id::Int)

    Worker(id, simulation.density, simulation.p_velocity, simulation.x_values, simulation.sources)
end
