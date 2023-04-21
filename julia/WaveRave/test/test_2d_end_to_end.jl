"""
Tests for running the 2D simulation end to end.
"""

using Test
using Debugger

using WaveRave

include("fixtures.jl")


@testset "2d end2end" begin
    @enter control = get_2d_test_control()
    end



