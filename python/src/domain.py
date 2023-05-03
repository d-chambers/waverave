#!/usr/bin/python3

import numpy as np
import os
import params
from . import simulation
from mpi4py import MPI
#import h5py




def make_homogeneous(nx, ny, dx, dy, vel):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
 
    par2d = simulation.simul_par(nx, ny , dx, dy, comm, rank, size)
    
    v = np.zeros((par2d.chunk_size, par2d.ny))
    v[:,:] = vel

    np.save('./outputs/vel_'+str(rank).zfill(5),v)


def make_layered(nx, ny, dx, dy, z, vel):

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
 
    par2d = simulation.simul_par(nx, ny , dx, dy, comm, rank, size)
    
    v = np.zeros((par2d.chunk_size, par2d.ny))
    for i,dep in enumerate(z):
        v[:dy*int(dep/dy)] = vel[i]

    np.save('./outputs/vel_'+str(rank).zfill(5),v)
    

