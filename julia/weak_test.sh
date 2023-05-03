#!/bin/bash
#SBATCH -N 12
#SBATCH -t 00:50:00
#SBATCH -j JWaveRaveWeak

# Record the node that we ran on
echo "Job ran on: $SLURM_NODELIST"

# Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1

# Load conda and activate julia environment
module load apps/python3/2022.10 
conda activate julia

# Run weaks scaling tests, adding workers increasing domain size in x direction
mpiexec -n 1 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,500
mpiexec -n 2 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,1000
# mpiexec -n 3 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,1500
# mpiexec -n 4 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,2000
# mpiexec -n 5 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,2500
# mpiexec -n 6 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,3000
# mpiexec -n 7 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,3500
# mpiexec -n 8 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,4000
# mpiexec -n 9 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,4500
# mpiexec -n 10 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,5000
# mpiexec -n 11 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,5500
# mpiexec -n 12 julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=300,6000
