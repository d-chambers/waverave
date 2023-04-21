"""
Module for simulations.
"""

using MPI

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
    MPI.Init()
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    size = MPI.Comm_size(MPI.COMM_WORLD)
    setup_io(rank, output_directory)
    local_sim = get_local_simulation(wave_sim, rank, size)




end