"""
A module for various utilities.
"""

using Combinatorics
using Memoization

export find_division_shape



"""
    Find the optimal shape for a given number of divisions.

Should work for any number of dimensions, but gets very expesnive for
dimensions above 3.
"""
@memoize function find_division_shape(divisions:: Int, domain_shape:: Tuple)
    
    shape_array = collect(domain_shape)

    # bail out if one division is requested
    if divisions == 1
        return ones(Int, length(shape_array)), shape_array
    end 
    
    # first get all possible shapes based on permutations that still have
    # a total of divisions
    ndims = length(domain_shape)
    perms = [digits(x - 1, base=divisions, pad=ndims) .+ 1 for x in 1:divisions^ndims]
    reduced_perms = Array{Int}(reduce(hcat, collect(perms))')
    valid_dims = [
        x for x in eachrow(reduced_perms)
        if (prod(x) == divisions) & (all((shape_array .% x) .== 0))
    ]
    if length(valid_dims) < 1
        msg = "$domain_shape cannot be divided into $divisions divisions"
        error(msg)
    end
    out = Array{Int}(reduce(hcat, valid_dims)')
    # now return the value that minmizes the sum of the dimensions
    chunks = shape_array' .÷ out
    (_, argmax) = findmin(sum(out, dims=2))
    return (out[argmax[1], :], chunks[argmax[1], :])
end


"""
    Chunk an array into n subdomains.

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
function chunk_grid(
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


function mkdir_if_not_exists(path:: String)
    if !isdir(path)
        mkdir(path)
    end
end


"""
    Function to determine if a location is inside the coordinates.
"""
function in_coords(location, coords)
    @assert length(location) == length(coords)
    in_coords = [
        (loc >= coor[1]) & (loc <= coor[end])
        for (loc, coor) in zip(location, coords)
    ]
    return all(in_coords)
end



"""
    Replace an element of an array, return a copy.
"""
function replace_ind(array, index::Int, value)
    out = copy(array)
    out[index] = value
    return out
end



"""
    Make homogeneous model
"""
function make_homogeneous_model(extents, dx, velocity)
    shape = [convert(Int, x ÷ dx) for x in extents]
    vel = convert(Float64, velocity)
    out = fill(vel, shape...)
    return out
end