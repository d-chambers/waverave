# Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1

# Load conda and activate julia environment
module load apps/python3/2022.10 
conda env create -n julia -c conda-forge -y

# activate julia environment and install all needed julia packages
conda activate julia
julia -e 'using Pkg; Pkg.activate("WaveRave"); Pkg.instantiate("WaveRave"); Pkg.resolve()'
