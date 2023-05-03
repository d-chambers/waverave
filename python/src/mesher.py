#!/usr/bin/python3

import sys
import numpy as np
import os
import params
import simulation
from mpi4py import MPI
#import h5py




def make_homogeneous(nx, ny, dx, dy, vel):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    par2d = simulation.simul_par(nx, ny , dx, dy, comm, rank, size)
    par2d.size_allocate()
    v = np.zeros((par2d.chunk_size, par2d.ny))
    v[:,:] = vel

    np.save(os.path.dirname(os.path.dirname(__file__))+'./outputs/vel_'+str(rank).zfill(5),v)


def make_layered(nx, ny, dx, dy, z, vel):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    print(rank)
    par2d = simulation.simul_par(nx, ny , dx, dy, comm, rank, size)
    par2d.size_allocate()
    

    v = np.zeros((par2d.chunk_size, par2d.ny))
    for i in range(1,len(z)):
        v[dy*int(z[i-1]/dy):dy*int(z[i]/dy)] = vel[i-1]
    #print(v)
    np.save(os.path.dirname(os.path.dirname(__file__))+'./outputs/vel_'+str(rank).zfill(5),v)
    return 0

if __name__ == "__main__":
    args = sys.argv
    nx  = int(args[1])
    ny  = int(args[2])
    dx  = params.dt
    dy  = params.dy
    z   = params.z
    vel = params.vel

    make_layered(nx, ny, dx, dy, z, vel)

   
