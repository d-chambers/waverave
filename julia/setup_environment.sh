# Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1

# Load conda and activate julia environment
module load apps/python3/2022.10 
conda create --yes --name julia --channel=conda-forge julia

# activate julia environment and install all needed julia packages
conda activate julia
julia -e 'using Pkg; Pkg.activate("WaveRave"); Pkg.instantiate(); Pkg.resolve()'
