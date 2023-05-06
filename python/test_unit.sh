#!/bin/bash
#SBATCH -N 2
#SBATCH --ntasks-per-node 6 
#SBATCH -t 00:15:00



module load mpi/openmpi/gcc/4.1.1
module load apps/python3/2022.10


conda activate mpirun

#mpirun -n 5 python test_unit.py 100 100

mpirun -n 5 python test_unit.py 100 100
