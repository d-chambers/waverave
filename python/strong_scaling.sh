#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node 12
#SBATCH -t 00:50:00
#SBATCH -J Pythonweak

# Record the node that we ran on
echo "Job ran on: $SLURM_NODELIST"

# Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1

# Load conda and activate mpirun environment
module load apps/python3/2022.10 
conda activate mpirun


#Run WeakScaling Tests

for i in {1..12}
do
mpirun -n $i python src/mesher.py 480 800
mpirun -n $i python src/run.py 480 800 5
done
