"""
Module defining sources.
"""


export AbstractSource, ExplosiveRickerSource, get_source_time_function, get_effective_frequencies

"""
Abstract type for a source. i.e. a Ricker wavelet.
"""
abstract type AbstractSource end


"""
A Ricker wavelet source from a shot.
    Source information.
    position should be a nd-array with the correct dimensions.
    time_delay controlls when the source fires in simulation time.
    frequency is in Hz.
    amplitude is the amplitude of the source.
"""
Base.@kwdef struct ExplosiveRickerSource <: AbstractSource
    location::AbstractArray
    time_delay::Real = 0.0
    frequency::Real = 20.0
    amplitude::Real = 1.0
end



"""
Get the source time function for a RickerSource.
"""
function get_source_time_function(source::ExplosiveRickerSource, time::AbstractArray{Real})
    τ_0 = time .- source.time_delay
    f_0 = source.frequency
    term1 = 1 - 2 * pi^2 * f_0^2 * τ_0^2
    term2 = exp(-pi^2 * f_0^2 * τ_0^2)
    return term1 * term2 * source.amplitude
end



"""
Get the "effective" frequencies for a source.
    We use 2.5 f_0 as the upper bound for the frequency, 
    which is consistent with specfem. We don't use the
    low-frequency yet so it is set as NaN
This is used to the wavelength 
"""
function get_effective_frequencies(source::ExplosiveRickerSource)::Tuple{Real,Real}
    f_0 = source.frequency
    return NaN, f_0 * 2.5
end
