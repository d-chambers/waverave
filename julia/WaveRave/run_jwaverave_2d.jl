"""
Script to run the Julia code for the WaveRave project using 2D inputs.
"""

using Base
using ArgParse


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "velocity_path"
            help = "The path to the velocity file (csv)"
            default = "velocity.csv"
        "source_path"
            help = "The path to the source file (csv)"
            default = "sources.csv"
        "--shape"
            help = "The shape to use for domain decomposition"
        "--out_path"
            help = "an option without argument, i.e. a flag"
            default = "out"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("$arg  => $val")
    end
end

main()
