"""
Script to run the Julia code for the WaveRave project using 2D inputs.
"""

using Base
using ArgParse
using WaveRave
using Debugger

include("test/fixtures.jl")


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "dx"
            help = "model spacing in m"
            default = 5
        "--extents"
            help = "The extents of the velocity model (z, x) in m"
            default = (300, 500)
        "--velocity"
            help = "The velocity for a homogeneous model (m/s)"
            default = 3000
        "--out_path"
            help = "where to save results."
            default = "outputs/results.hdf5"
        "--source_frequency"
            help = "frequency of Ricker wavelet in Hz"
            default = 20.0
        "--source_delay"
            help = "Source delay in s"
            default = 0.5
        "--save_interval"
            help = "interval to save the wavefield in s"
            default = 1.0
        "--simulation_time"
            help = "total time (in s) to run simulation"
            default = 5.0
        "--runs"
            help = "total number of runs for averaging times"
            default = 4
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    # ensure expected datatypes from each command line argument
    extents = [convert(Float64, x) for x in parsed_args["extents"]]
    dx = convert(Float64, parsed_args["dx"])
    vel = convert(Float64, parsed_args["velocity"])
    source_freq = convert(Float64, parsed_args["source_frequency"])
    source_delay = convert(Float64, parsed_args["source_delay"])
    save_interval = convert(Float64, parsed_args["save_interval"])
    sim_time = convert(Float64, parsed_args["simulation_time"])
    # init velocity model
    velmod = make_homogeneous_model(extents , dx, vel)
    coords = [collect(0:dx:(x-dx)) .+ dx/2 for x in extents]
    # init source
    source = ExplosiveRickerSource(
        location=[x/2 for x in extents],
        frequency=source_freq,
        time_delay=source_delay,
    )
    # init wavesim
    sim = WaveSimulation(
        coords=coords,
        p_velocity=velmod,
        sources=[source],
        time_max=sim_time,
    )
    sim()  # validate simulation params
    run_wave_simulation(
        sim; 
        snapshot_interval=save_interval,
        save_path=parsed_args["out_path"],
    )
end

main()
