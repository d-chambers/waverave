module WaveRave

using Reexport

include("Stencil.jl")
export get_stencil, get_stencil_coefs

include("Control.jl")
export WaveSimulation, validate_simulation, LocalGrid

include("Sources.jl")
export ExplosiveRickerSource, get_source_time_function

include("Utils.jl")
export find_division_shape, make_homogeneous_model

include("Decompose.jl")
export get_domain_map

include("Simulation.jl")
export run_wave_simulation

end
