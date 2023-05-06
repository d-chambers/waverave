## Import modules
module load mpi/openmpi/gcc/4.1.1
module load apps/python3/2022.10


##create conda environment

conda create --name mpirun python=3.8 numpy matplotlib mpi4py


