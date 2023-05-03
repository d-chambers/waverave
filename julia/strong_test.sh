#!/bin/bash
#SBATCH -N 12
#SBATCH -t 00:50:00
#SBATCH -J JWRStrong

# Record the node that we ran on
echo "Job ran on: $SLURM_NODELIST"

# Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1

# Load conda and activate julia environment
module load apps/python3/2022.10 
conda activate julia

# Run julia script
for n in {1..12} 
do
    mpiexec -n $n julia --project=WaveRave WaveRave/run_jwaverave_2d.jl --profile --extents=2400,4000
done
