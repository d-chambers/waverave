"""
Functions for decomposing domains into subdomains.
"""

using Base
using Combinatorics


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
