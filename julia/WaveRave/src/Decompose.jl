"""
Functions for decomposing domains into subdomains.
"""

using Base
using Combinatorics
using ImageFiltering

include("Control.jl")
include("Utils.jl")


"""
    Decompose a 1D array into n subdomains.

# Arguments
- `array:: AbstractArray{Real, 1}`: The array to be decomposed.
- `x_values:: AbstractArray{Real, 1}`: The x values for the array.
- `id:: Int`: The ID of the subdomain to be returned.
- `id_max:: Int`: The total number of subdomains.
- `overlap:: Int`: The number of points of overlap between subdomains.
- `boundary_fill:: Symbol`: The method for filling the boundary of the 
    subdomain. Options are:
    - `:wraparound`: The boundary is filled with the opposite side of the 
        array.
    - `:zero`: The boundary is filled with zeros.
    - `:repeat`: The boundary is filled with the nearest value.
"""
function grid_decomposition(
    array:: AbstractArray{Real, 1}, 
    x_values:: AbstractArray{AbstractArray{Real, 1}}, 
    id:: Int,
    id_max:: Int, 
    overlap:: Int,
    boundary_fill:: Symbol = :wraparound
    ):: Tuple(AbstractArray{Real, 1}, AbstractArray{Real, 1})

    _check_inputs(array, x_values, id, id_max, overlap, boundary_fill)
    
end


function _check_inputs(
    array, 
    x_values, 
    id:: Int,
    id_max:: Int, 
    overlap:: Int,
    boundary_fill:: Symbol,
)
    # basic input checks
    @assert id <= id_max "id must be less than or equal to id_max"
    @assert id > 0 "id must be greater than 0"  
    @assert id_max > 0 "id_max must be greater than 0"
    @assert overlap >= 0 "overlap must be greater than or equal to 0"
    @assert boundary_fill ∈ (:wraparound, :zero, :repeat) "Invalid fill method"
    # check that the array and x_values are the same length
    @assert size(array) == [len(x) for x in x_values] "array and x_values must be the same length"

end


"""
    Get the length of the division, before padding
"""
function get_div_len(rank, size, coord_len)
    div_len = coord_len ÷ size
    if coord_len % size > rank  # if there is leftover sizes
        div_len += 1
    end
    return div_len
end

using Debugger
"""
    Fill the array with closest values around padding.
"""
function fill_pad_array!(array, pad)
    # first get corners
    combs = [(1:pad + 1, x-(pad+1):x) for x in size(array)]

    for inds in Iterators.product(combs...)
        vals = array[inds...]
        array[inds...] .= maximum(vals[isfinite.(vals)] )
    end
    # then get rows/columns
    for i in pad:size(array)[1]-pad
        array[i, 1:pad] .= array[i, pad + 1] 
        array[i, end-pad:end] .= array[i, end - pad - 1]
    end
    for i in pad:size(array)[2]-pad
        array[1:pad, i] .= array[pad + 1, i] 
        array[end-pad:end, i] .= array[end - pad - 1, i]
    end

    return array
end


"""
    Decompose an array into n subdomains using pencil decomposition.

Returns both the new coords and the padded arrays.
Decomposes over the largest dimension and fills the padded regions
outside of the grid with the last grid value. Assumes 2D array
"""
function pencil_decomposition(coords, array, pad, rank, rank_count)
    @assert(pad >= 0, "pad must be greater than or equal to 0")
    # get the largest dimension
    @assert(length(coords) == 2, "only 2d for now.")
    origin_index = [1:length(x) for x in coords]
    dim = argmax(size(array))
    out_coords = copy(coords)
    div_coord = out_coords[dim]
    cord_len = length(div_coord)

    # get div lens    
    div_lens = [get_div_len(x, rank_count, cord_len) for x in 0:1:rank_count-1]
    interfaces = cat([1], cumsum(div_lens), dims=1)
    start_ind, end_ind = interfaces[rank + 1], interfaces[rank + 2]
    origin_index[dim] = start_ind:end_ind
    # get start of division
    out_coords[dim] = div_coord[start_ind: end_ind]
    # init padded array and fill with values from original array
    pad_arrays = [x.start + pad:pad + x.stop for x in origin_index]
    padded_size = [length(x) for x in out_coords] .+ 2 * pad
    out = padarray(array[origin_index...], Pad(:replicate,fill(pad, length(coords))...))
    # last assignment ensures nothing was overwritten with sloppy edge padding.
    # TODO: We need to make sure to update the values after padding with 
    # correct values from adjacent local blocks.
    return out_coords, out, origin_index
end





"""
Get the local simulation for a given rank of specified size.
"""
function get_local_simulation(
    global_simulation:: WaveSimulation,
    rank:: Int,
    rank_count:: Int,
    decomposition:: Symbol = :pencil,
    ) :: LocalGrid
    if decomposition != :pencil
        error("only pencil decomposition supported now.")
    end

    coords, vel, global_index = pencil_decomposition(
        global_simulation.coords, 
        global_simulation.p_velocity, 
        global_simulation.space_order,
        rank,
        rank_count,
    )

    sources = [
        x for x in global_simulation.sources 
        if in_coords(x.location, coords)
    ]
    receivers = [
        x for x in global_simulation.receivers 
        if in_coords(x.location, coords)
    ]

    out = LocalGrid(
        coords=coords,
        p_velocity=vel,
        sources=sources,
        receivers=receivers,
        global_index=global_index,
    )
    return out
end
    
