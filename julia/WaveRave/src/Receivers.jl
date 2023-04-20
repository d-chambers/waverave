"""
Module for handling receivers.
"""


"""
    Class to define seismic Receivers.
"""
Base.@kwdef struct Receiver
    name:: String
    location:: AbstractArray
end
