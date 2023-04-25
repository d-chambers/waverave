"""
Tests for running the 2D simulation end to end.
"""

using Test
using Debugger
using WaveRave

include("fixtures.jl")


@testset "2d end2end serial" begin
    global_sim = get_2d_test_control()()
    @enter run_wave_simulation(global_sim)

    end



