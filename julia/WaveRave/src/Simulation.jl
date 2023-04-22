"""
Module for simulations.
"""

using MPI
using Debugger

include("Control.jl")
include("Utils.jl")
include("Decompose.jl")

"""
Ensure the output directory exists for all ranks.
"""
function setup_io(rank, output_directory)
    if rank == 0
        mkdir_if_not_exists(output_directory)
    end
    MPI.Barrier(MPI.COMM_WORLD)
end


"""
    Run the simulation.
# Arguments
- wave_sim - Parameters of the simulation to run.
- snapshot_interval - The interval at which to save snapshots of the simulation (in seconds)
- output_directory - The directory to save the snapshots to.
"""
function run_wave_simulation(
    wave_sim:: WaveSimulation, 
    boundary_condition:: Symbol = :reflective,
    snapshot_interval:: Union{Int, Nothing} = nothing,
    output_directory:: String = "outputs",
    )
    # MPI setup
    MPI.Init()
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    rank_count = MPI.Comm_size(MPI.COMM_WORLD)
    # create output directory and get local version
    setup_io(rank, output_directory)
    local_sim::LocalGrid = get_local_simulation(wave_sim, rank, rank_count)
    dx = [x[2] - x[1] for x in local_sim.coords]  # dx in all dims
    coefs = get_stencil_coefs(wave_sim.space_order, derivative=2)
    # initiliaze arrays for previous, current, and next time.
    tp = zeros(shape(local_sim.p_velocity))
    tc = zeros(shape(local_sim.p_velocity))
    tn = zeros(shape(local_sim.p_velocity))
    
    if 
    
end