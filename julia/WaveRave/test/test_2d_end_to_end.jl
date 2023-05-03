"""
Tests for running the 2D simulation end to end.
"""

using Base
using Test
using Debugger
using JLD2
using MPI

using WaveRave


include("fixtures.jl")


# Tests to ensure MPI with different numbers of nodes returns the same wavefield
@testset "2d end2end mpi" begin
    jwaverave_path = dirname(Base.active_project())
    script_path = joinpath(jwaverave_path, "run_jwaverave_2d.jl")
    out_path = joinpath(joinpath(joinpath(jwaverave_path, "test"), "outputs"), "mpi_comps")
    julia_path = Base.julia_cmd()
    n = 2  # number of processess
    sims = Dict()
    mpiexec() do exe  # MPI wrapper
        for n in 1:4
            opath = joinpath(out_path, "mpi_$n.hdf5")
            println("Writing $opath in MPI test")
            cmds = Cmd(
                `mpiexec -n $n $julia_path --project=$jwaverave_path $script_path --savepath=$opath`,
                )
            out = run(cmds)
            @test isfile(opath)
            sims[n] = JLD2.load_object(opath)
        end
    end
    # read each of the created files and ensure they are the same.
    for (n, sim1) in sims
        for (n, sim2) in sims
            @assert sim1.wavefield ≈ sim2.wavefield
        end
    end
end



@testset "2d end2end serial" begin
    global_sim = get_2d_test_control()()
    path = "outputs/2d_serial_wavesim.hdf5"
    (sim, _, _) = run_wave_simulation(global_sim; save_path=path, snapshot_interval=0.1)
    # wavefield should be non_empty
    @test maximum(sim.wavefield) > 0
    # output file should have been created
    @test isfile(path)
    # saved file should round-trip
    sim2 = JLD2.load_object(path)
    @test sim.wavefield == sim2.wavefield
    @test sim.coords == sim2.coords
end
    

# end

