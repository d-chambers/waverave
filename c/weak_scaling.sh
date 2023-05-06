#!/bin/bash
#SBATCH -N 12
#SBATCH --ntasks-per-node 12
#SBATCH -t 00:50:00

# Record the node that we ran on
echo "Job ran on: $SLURM_NODELIST"

# Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1
mpicc mpi_simul.c -o mpi_simul -lm
for n in {1..12}
do
	mpirun -np $n ./mpi_simul  $n*100 60
done

