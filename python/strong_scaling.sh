#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node 12
#SBATCH -t 00:50:00
#SBATCH -J Pythonstrong

# Record the node that we ran on
echo "Job ran on: $SLURM_NODELIST"

# Load module for mpi (assuming gcc loaded)
module load mpi/openmpi/gcc/4.1.1

# Load conda and activate mpirun environment
module load apps/python3/2022.10 
conda activate mpirun


#Run WeakScaling Tests
mpirun -n 1 python src/mesher.py 480 800
mpirun -n 1 python src/run.py 480 800 5
mpirun -n 2 python src/mesher.py 480 800
mpirun -n 2 python src/run.py 480 800 5
mpirun -n 3 python src/mesher.py 480 800
mpirun -n 3 python src/run.py 480 800 5
mpirun -n 4 python src/mesher.py 480 800
mpirun -n 4 python src/run.py 480 800 5
mpirun -n 5 python src/mesher.py 480 800
mpirun -n 5 python src/run.py 480 800 5
mpirun -n 6 python src/mesher.py 480 800
mpirun -n 6 python src/run.py 480 800 5
mpirun -n 7 python src/mesher.py 480 800
mpirun -n 7 python src/run.py 480 800 5
mpirun -n 8 python src/mesher.py 480 800
mpirun -n 8 python src/run.py 480 800 5
mpirun -n 9 python src/mesher.py 480 800 
mpirun -n 9 python src/run.py 480 800 5
mpirun -n 10 python src/mesher.py 480 800 
mpirun -n 10 python src/run.py 480 800 5
mpirun -n 11 python src/mesher.py 480 800 
mpirun -n 11 python src/run.py 480 800 5
mpirun -n 12 python src/mesher.py 480 800
mpirun -n 12 python src/run.py 480 800 5
