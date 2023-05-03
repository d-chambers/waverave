"""
Script to make a plot of the output wavefields
"""

using WaveRave
using JLD2
using Debugger
using Printf
using Plots

global data_path = "outputs/results.hdf5"
global plot_path = "outputs/wavefield_plots"


function mkdir_if_not_exists(path:: String)
    if !isdir(path)
        mkdir(path)
    end
end


function plot_wavefield_time_slice(wf, coords, time; clim=nothing)
    (zs, xs) = coords
    heatmap(xs, zs, wf, c=cgrad([:blue,:white,:red]), clim=clim, legend=:none)

    time_str = @sprintf("%.03f", time)
    plot!(title = "time=$time_str", xlabel = "x (m)", ylabel = "z (m)")
    # plot!(xlims=xlims, ylims=zlims)
    save_path = "$plot_path/$time_str.png"
    savefig(save_path)

    
end


function main()
    max_scalar = 4
    obj = JLD2.load_object(data_path)
    rm(plot_path, force=true, recursive=true)
    mkdir_if_not_exists(plot_path)
    time_vector = obj.time_vector
    wf = obj.wavefield
    abs_max = maximum(abs.(wf))
    clim = (-abs_max/max_scalar, abs_max/max_scalar)
    for (ind, time) in enumerate(time_vector)
        sub = wf[ind, :, :]
        plot_wavefield_time_slice(sub, obj.coords, time; clim=clim)

    end
    
end

main()