module WaveRave
    include("Stencil.jl")
    export get_stencil

    include("Control.jl")
    export WaveSimulation, validate_simulation

    include("Sources.jl")
    export ExplosiveRickerSource, get_source_time_function
end
