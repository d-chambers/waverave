"""
A module for various utilities.
"""

using Combinatorics

export find_division_shape

function atleast2d(x::AbstractArray{<:Any, 1})
    return reshape(x, length(x), 1)
end


function atleast2d(x::Any)
    return x
end




"""
    Find the optimal shape for a given number of divisions.

Should work for any number of dimensions, but gets very expesnive for
dimensions above 3.
"""
function find_division_shape(divisions:: Int, domain_shape:: Tuple)
    
    shape_array = collect(domain_shape)

    # bail out if one division is requested
    if divisions == 1
        return ones(Int, length(shape_array)), shape_array
    end 
    
    # first get all possible shapes based on permutations that still have
    # a total of divisions
    ndims = length(domain_shape)
    perms = [digits(x - 1, base=divisions, pad=ndims) .+ 1 for x in 1:divisions^ndims]
    reduced_perms = reduce(hcat, collect(perms))'
    sub_1 = reduced_perms[prod(reduced_perms, dims=2) .== divisions]
    sub_1_2d = atleast2d(sub_1)
    # next, eliminat any solutuions that don't have a int number of elements
    sub_2 = filter(x -> all((shape_array .% x) .== 0), sub_1_2d)
    sub_2_2d = atleast2d(sub_2)
    if length(sub_2) < 1
        msg = "$domain_shape cannot be divided into $divisions divisions"
        error(msg)
    end
    # now return the value that maximizes block volume (minimizes ghost zones)
    (_, argmax) = findmax(prod(sub_2_2d, dims=2))
    best_divisor = sub_2_2d[argmax[1], :]
    return (best_divisor, shape_array .÷ best_divisor)
    
end

