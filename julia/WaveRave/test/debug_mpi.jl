"""
A script to debug MPI stuff.
"""

using WaveRave
using Debugger

function main()
    dx = 5
    extents = (300, 500)
    vel = 3000
    source_freq = 60.0
    source_delay = 0.1
    save_interval = 0.1
    sim_time = 0.5
    space_order=4



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
        space_order=space_order,
    )
    sim()  # validate simulation params

    @enter run_wave_simulation(
        sim; 
        snapshot_interval=save_interval,
        _fake_rank = 1,
        _fake_rank_count = 4,
    )
end

main()