"""
Module for constants used in the WaveRave package.
"""


"""
A type for a 1,2, or 3 dimensional array of Real numbers.
"""
OneToThreeDimRealArray = Union{
    AbstractArray{Real, 3}, 
    AbstractArray{Real, 2}, 
    AbstractArray{Real, 1},
}
