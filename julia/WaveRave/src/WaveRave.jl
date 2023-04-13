module WaveRave
    include("Stencil.jl")
    include("Control.jl")

    export get_stencil, WaveSimulation, RickerSource
end
