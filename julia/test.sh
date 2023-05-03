# run the Julia package tests for WaveRave
julia --project="./WaveRave" -e 'include("WaveRave/test/runtests.jl")'
# julia --project="./WaveRave" -e "using Base; println(Base.active_project())"
