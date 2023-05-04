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
   
    par2d = simulation.simul_par(nx, ny , dx, dy, comm, rank, size)
    par2d.size_allocate()

    v = np.zeros((par2d.chunk_size, par2d.ny))
    v[:,:] = vel[-1]
    for i in range(2, len(z)):
        if  z[-i] < par2d.end and z[-i] >= par2d.start:
            if z[-i-1]>= par2d.start and z[-i-1] < par2d.end:
               start = z[-i-1] - par2d.start
            else:
                start = 0
            end = z[-i] - par2d.start
            v[start:end,:] = vel[-i]
        elif par2d.end <= z[-i]:
             v[:,:] = vel[-i]
           


                 
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

    make_homogeneous(nx, ny, dx, dy, vel)

   
