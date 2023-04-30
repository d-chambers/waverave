#!/bin/bash
#SBATCH -N 3 
#SBATCH --ntasks-per-node=4
#SBATCH -t 00:01:00

# Record the node that we ran on
echo "Job ran on: $SLURM_NODELIST"

# Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1

# Load conda and activate julia environment
module load apps/python3/2022.10 
conda activate julia

# Run julia script
mpiexec -n 4 julia mpi_test.jl
