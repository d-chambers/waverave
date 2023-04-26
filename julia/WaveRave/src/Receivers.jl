"""
Module for handling receivers.
"""

module Receivers

export Receiver

"""
    Class to define seismic Receivers.
"""
Base.@kwdef struct Receiver
    name:: String
    location:: AbstractArray
end

end