"""
Module for simulations.
"""

using MPI
using Debugger
using ImageFiltering

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
    Get the index in a local grid. Must not be in padded region.

Simply finds the closest point, but can't be far from the edge.
"""
function get_inds(ar, coords)
    out = [
        argmin(abs.(coord .- loc)) 
        for (coord, loc) in zip(coords, ar)
    ]
    return out

end


"""
Exchange data between processes 
"""
function update_ghost_regions(array, domain_map, rank)
    rank_count = length(domain_map.local_coord_map)
    if rank_count == 1  # only on process, do nothing
        return
    end

end


"""
Make Laplacian image
"""
function make_laplace_image(wave_sim)
    image = zeros(fill(wave_sim.space_order*2+1, length(wave_sim.coords))...)
    center_coords::Array{Any} = fill(wave_sim.space_order+1, length(wave_sim.coords)...)
    for (num, coord) in enumerate(wave_sim.coords)
        imge_coords = copy(center_coords)
        imge_coords[num] = Colon()
        dx = coord[2] - coord[1]
        stencil = get_stencil_coefs(wave_sim.space_order) ./ (dx^2)
        image[imge_coords...] .+= stencil
    end
    return image
    
end


"""
    Apply boundary conditions
"""
function apply_boundary(array, domain_map, rank, type=:reflective)
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
    # basic setup
    time = collect(0:wave_sim.dt:wave_sim.time_max)
    setup_io(rank, output_directory)
    domain_map::DomainMap = get_domain_map(wave_sim, rank_count)
    # get padded p and update with values from other ranks
    local_p = get_padded_local_array(domain_map, wave_sim.p_velocity, rank)
    update_ghost_regions(local_p, domain_map, rank)
    local_sources, local_receivers = get_local_sources(wave_sim, domain_map, rank) 
    coords = domain_map.local_coord_map[rank]
    dt_sq = wave_sim.dt ^ 2
    vel_sq = local_p .^ 2
    image = make_laplace_image(wave_sim)
    # initiliaze arrays for previous, current, and next time.
    tp = zeros(size(local_p))
    tc = zeros(size(local_p))
    # get the index for the closest grid point to source.
    source_inds_no_pad = [get_inds(x.location, coords) for x in local_sources]
    source_inds = [x .+ wave_sim.space_order for x in source_inds_no_pad]
    source_time_functions = [x(time) for x in local_sources]
    # time step
    for (t_ind, t) in enumerate(time)
        update_ghost_regions(tc, domain_map, rank)
        # save 
        if !isa(snapshot_interval, Nothing) && (t_ind % snapshot_interval == 0)
            save_wave(tc)
        end
        @bp
        # run simulation
        laplace = imfilter(tc, image) # , Inner())
        tn = @. vel_sq * laplace * dt_sq + 2 * tc - tp
        for (source_ind, stf) in zip(source_inds, source_time_functions)
            tn[source_ind...] .+= stf[t_ind] .* dt_sq 
        end
        apply_boundary(tn, domain_map, rank, boundary_condition)
        # rotate time arrays
        tc, tp = tn, tc
    end
    
end