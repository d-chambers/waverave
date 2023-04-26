"""
Module for getting stencils of varius lengths
"""

module Stencil

export get_stencil, get_stencil_coefs

"""
    get_stencil(x::AbstractVector{<:Real}, x₀::Real, m::Integer)

A function to get the stencil of length `m` for the point `x₀` in the
grid `x`. The stencil is a vector of weights `w` such that the `m`the
derivative of a function `f` at `x₀` is approximated by `w ⋅ f(x)`.

The code was shamelessly stolen from 
https://discourse.julialang.org/t/generating-finite-difference-stencils/85876/5
"""
function get_stencil(x::AbstractVector{<:Real}, x₀::Real, m::Integer)
    ℓ = 0:length(x)-1
    m in ℓ || throw(ArgumentError("order $m ∉ $ℓ"))
    A = @. (x' - x₀)^ℓ / factorial(ℓ)
    return A \ (ℓ .== m) # vector of weights w
end


"""
    Get generic stencil coeficients for central difference.
"""
function get_stencil_coefs(order, derivative=2)
    dx = -order:1:order
    return get_stencil(dx, 0, derivative)
end

end