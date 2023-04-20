"""
Tests for the control module.
"""

using WaveRave
using Test
using Debugger

include("fixtures.jl")


# ---- Tests -----


@testset "test 1d" begin
    @testset "simple init" begin
        wc = get_simple_1d_control()
        @test isa(wc, WaveSimulation)
    end

    @testset "test CFL limit" begin
        # case that should pass
        wc = get_simple_1d_control()
        @test validate_simulation(wc)
        # case that should fail
        wc = get_simple_1d_bad_clf()
        @test_throws "CFL limit exceeded!" validate_simulation(wc)
    end

    @testset "test calling control validates it." begin
        wc = get_simple_1d_control()()
        @test isa(wc, WaveSimulation)
    end

end
