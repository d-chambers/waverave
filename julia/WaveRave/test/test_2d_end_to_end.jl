"""
Tests for running the 2D simulation end to end.
"""

using Test
using Debugger
using JLD2
using MPI

using WaveRave


include("fixtures.jl")

function see(sim1, sim2)
    println(sim)
end

@testset "2d end2end serial" begin
    global_sim = get_2d_test_control()()
    path = "outputs/2d_serial_wavesim.hdf5"
    sim = run_wave_simulation(global_sim; save_path=path, snapshot_interval=0.1)
    # output file should have been created
    @test isfile(path)
    # saved file should round-trip
    sim2 = JLD2.load_object(path)
    @test sim.wavefield == sim2.wavefield
    @test sim.time_vector == sim2.time_vector
    @test sim.coords == sim2.coords
end


@testset "2d end2end mpi" begin
    # global_sim = get_2d_test_control()()
    # path = "outputs/2d_2proc.hdf5"
    # n = 2  # number of processes
    # mpiexec() do exe  # MPI wrapper
    #     run(`$exe -n $n $(Base.julia_cmd()) run_jwaverave.jl --out_path=$path`)
    # end

end

