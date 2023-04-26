"""
Module for simulations.
"""
module Simulation

using ..Decompose
using ..Control
using ..Utils
using MPI
using Debugger
using ImageFiltering
import JLD2


export run_wave_simulation


"""
Ensure the output directory exists for all ranks.
"""
function setup_io(rank, save_path)
    if rank == 0
        path = isdir(save_path) ? save_path : dirname(save_path)
        mkdir_if_not_exists(path)
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
    if rank_count == 1  # only one process, do nothing
        return
    end
    neighbors = domain_map[rank]
end


"""
    Put the global wavefields back together.
"""
function get_global_wavefield(wave_field, domain_map::DomainMap, rank)
    # setup rank 0 to receive over wavefields
    rank_count = length(domain_map.local_coord_map)
    ranks = collect(0:(rank_count -1))
    comm = MPI.COMM_WORLD
    # Fast exit ramp for single process
    if rank_count == 1  
        return wave_field
    end
    
    requests::AbstractArray{MPI.Request} = []
    if rank == 0
        rec_dict:: Dict{Int=>AbstractArray} = Dict()
        glob_inds = [x[1] for x in domain_map.global_inds]
        glob_size = vcat(size(wave_feild)[1:1], glob_inds)
        out = zeros(Float64, glob_size...)
        for r in ranks
            recv_mesg = Array{Float64}(undef, ndims(out))
            push!(requests, MPI.Irecv!(recv_mesg, comm; source=r, tag=r+32))
            rec_dict[r] = recv_mesg
        end
        # wait for all messages to be received by rank 0
        MPI.Waitall(requests)
        for (r, data) in rec_dict
            inds = domain_map.local_index_map[r]
            out[inds...] = data
        end
        return out
    else
        push!(requests, MPI.Isend(wave_field, comm; dest=0, tag=rank+32))
        return []
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
    Apply boundary conditions.

This works by first finding the domain edges that arent shared, then
applying the 0 condition.
"""
function apply_boundary(array, domain_map::DomainMap, rank, type=:reflective)
    @assert type == :reflective
    nmap = domain_map.neighbors_map[rank]
    ndim = length(domain_map.global_coords)
    inds = domain_map.local_index_map[rank]
    pad = domain_map.padding

    colons::AbstractArray{Any} = [Colon() for _ in 1:ndim]
    start_inds = [1:1+pad for _ in inds]
    ind_max = [x[2] + pad: x[2] + 2*pad for x in eachrow(inds)]
    
    for (dim, row) in enumerate(eachrow(nmap))
        for (ind_range, val) in zip([start_inds, ind_max], row)
            # this dimension has a neighbor; keep going.
            if ~isa(val, Nothing)
                continue
            end
            # otherwise set values to first non-padded index to 0
            new_ind = replace_ind(colons, dim, ind_range[dim])
            array[new_ind...] .= 0
        end
    end
    

end


"""
    Init the output wavefield.
"""
function init_output_wavefield(time, domain_map::DomainMap, rank, snapshot)
    tmax = maximum(time)
    potential_saves = collect(0:snapshot:(tmax ÷ snapshot)*snapshot)
    save_times = [time[searchsortedfirst(time, x)] for x in potential_saves]
    inds = domain_map.local_index_map[rank]
    wave_inds = vcat([length(save_times)], [x[2] for x in eachrow(inds)])
    out = zeros(Float64, wave_inds...)
    return out, Set(save_times)
end


"""
    Run the simulation.
# Arguments
- wave_sim - Parameters of the simulation to run.
- snapshot_interval - The interval at which to save snapshots of the simulation (in seconds)
- output_directory - The directory to save the snapshots to.
"""
function run_wave_simulation(
    wave_sim:: WaveSimulation; 
    boundary_condition:: Symbol = :reflective,
    snapshot_interval::Float64 = 1.0,
    save_path:: String = "outputs/wavesim.hdf5",
    ) :: WaveSimulation
    # MPI setup
    MPI.Init()
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    rank_count = MPI.Comm_size(MPI.COMM_WORLD)
    # basic setup
    time = collect(0:wave_sim.dt:wave_sim.time_max)
    setup_io(rank, save_path)
    domain_map::DomainMap = get_domain_map(wave_sim, rank_count)
    (wavefield, save_times) = init_output_wavefield(
        time, domain_map, rank, snapshot_interval
    )
    # get padded p and update with values from other ranks
    local_p = get_padded_local_array(domain_map, wave_sim.p_velocity, rank)
    local_unpad_coords = domain_map.unpad_coord_map[rank]
    update_ghost_regions(local_p, domain_map, rank)
    local_sources, local_receivers = get_local_sources(wave_sim, domain_map, rank) 
    coords = domain_map.local_coord_map[rank]
    dt_sq = wave_sim.dt ^ 2
    vel_sq = local_p .^ 2
    image = centered(make_laplace_image(wave_sim))
    # initiliaze arrays for previous, current, and next time.
    tp = zeros(size(local_p))
    tc = zeros(size(local_p))
    # get the index for the closest grid point to source.
    source_inds_no_pad = [get_inds(x.location, coords) for x in local_sources]
    source_inds = [x .+ wave_sim.space_order for x in source_inds_no_pad]
    source_time_functions = [x(time) for x in local_sources]
    update_count = 1
    # time step
    for (t_ind, t) in enumerate(time)
        update_ghost_regions(tc, domain_map, rank)
        # save 
        if t ∈ save_times
            new_ind = vcat([update_count], [Colon() for _ in 1:ndims(tc)])
            wavefield[new_ind...] = tc[local_unpad_coords...]
            update_count += 1
        end
        # run simulation
        laplace = imfilter(tc, image) # , Inner())
        tn = @. vel_sq * laplace * dt_sq + 2 * tc - tp
        # add source contributions
        for (source_ind, stf) in zip(source_inds, source_time_functions)
            tn[source_ind...] += (stf[t_ind] * dt_sq) 
        end
        apply_boundary(tn, domain_map, rank, boundary_condition)
        # rotate time arrays
        tc, tp = tn, tc
    end
    # get global wavefield and update output.
    global_wavefield = get_global_wavefield(wavefield, domain_map, rank)
    wave_sim.wavefield = global_wavefield
    wave_sim.time_vector = sort(collect(save_times))
    if rank == 0
        JLD2.save_object(save_path, wave_sim)
    end
    return wave_sim
end

end
