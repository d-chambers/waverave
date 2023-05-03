"""
Tests for decomposing domains.
"""

using Test
using Debugger
using WaveRave

include("fixtures.jl")

function something(dm)
    a = 1
end

@testset "1D grid_decomposition" begin
    @testset "1D 2" begin
        # Test 1D decomposition.
        control = get_simple_1d_control()
        dm = get_domain_map(control, 2)
        # there should be two domains
        @test length(dm.local_coord_map) == 2
        @test length(dm.local_index_map) == 2
        @test length(dm.local_coord_limit_map) == 2        
    end
    @testset "1D 3" begin
        # Test 1D decomposition.
        control = get_simple_1d_control()
        dm = get_domain_map(control, 3)
        # there should be two domains
        @test length(dm.local_coord_map) == 3
        @test length(dm.local_index_map) == 3
        @test length(dm.local_coord_limit_map) == 3
    end
end


@testset "2D grid_decomposition" begin
    @testset "4 divisions" begin
        # Test 2D decomposition.
        control = get_2d_test_control()
        dm = get_domain_map(control, 4)
        # there should be two domains
        @test length(dm.local_coord_map) == 4
        @test length(dm.local_index_map) == 4
        @test length(dm.local_coord_limit_map) == 4
    end
    @testset "2 divisions" begin
        # Test 2D decomposition.
        control = get_2d_test_control()
        @enter get_domain_map(control, 2)
        dm = get_domain_map(control, 2)
        # there should be two domains
        @test length(dm.local_coord_map) == 4
        @test length(dm.local_index_map) == 4
        @test length(dm.local_coord_limit_map) == 4
    end

end