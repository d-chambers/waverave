#!/usr/bin/python3

import sys
import solver
import simulation 
import numpy as np
from mpi4py import MPI
import params
import os
import time


def time_simulation(start, end, wavefield, rank_source, isx, isy, F,  dt, output_dir, comm):
    """
    Run Time Simulations for given time steps
    start: Starting time step
    end:   Ending Time Step
    wavefield: A wavefield class object
    rank_source: Rank of the wavefield containing source point
    isx: Source Location in x direction
    isy: Source Location in y direction
    dt: Time step 
    """
    for i in range(start, end):
        if wavefield.rank == rank_source:
            wavefield.inject_source(isx, isy, F[i], dt)
        recv_status = []
        send_status = []
        recv_request, send_request = wavefield.ghost_region()

        MPI.Request.waitall(recv_request, recv_status)
        MPI.Request.waitall(send_request, send_status)
        comm.Barrier()



        wavefield.update_solution()
        wavefield.set_boundaries()

        tmp = wavefield.new
        wavefield.new = wavefield.old
        wavefield.old = tmp

        
        
    wavefield.write(end, output_dir)


def time_simulation_with_checkpoint(start, stop, wavefield, rank_source, isx, isy, F, dt, checkpoint, output_dir, comm):
    """
    Time solution with checkpointing
    Run Time Simulations for given time steps
    start: Starting time step
    end:   Ending Time Step
    wavefield: A wavefield class object
    rank_source: Rank of the wavefield containing source point
    isx: Source Location in x direction
    isy: Source Location in y direction
    dt: Time step
    checkpoint: Checkpoint Steps
    """
    time = np.arange(start, stop, checkpoint)
    if stop not in time:
        time = np.append(time,stop)
    for i in range(len(time)-1):
        time_simulation(time[i], time[i+1], wavefield, rank_source, isx, isy, F, dt, output_dir,comm)



def solve_2d_mpi(nx, ny, nt, sx, sy, dx, dy, dt, pad, type, checkpoint, output_dir):
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if type == 0:
       source = solver.Ricker(ricker_sigma, ricker_shift, nt, dt)
       F = source.source_time()
    else:
         print("No other source implemented other than Ricker")
         source = solver.Ricker(nt, dt)
         F = source.source_time()


    par2d = simulation.simul_par(nx, ny , dx, dy, comm, rank, size)
    par2d.size_allocate()
    rank_source, isx, isy = par2d.source(sx, sy, pad)
    
    

    file = os.path.dirname(os.path.dirname(__file__))+'/outputs/vel_'+str(rank).zfill(5)+'.npy'

    wavefield = solver.wavefield(par2d.lnx, par2d.lny, pad, dx, dy, dt, comm, rank, size, par2d.method)
    wavefield.vel(file)
    wavefield.set_exchange()

    start=time.time()
    if checkpoint == 0:
       time_simulation(0, nt, wavefield, rank_source, isx, isy, F ,  dt, output_dir, comm) 
       
    
    else:
         time_simulation_with_checkpoint(0, nt, wavefield, rank_source, isx, isy, F, dt, checkpoint, output_dir, comm)
    if rank==0:
       end=time.time()
       print("Average Time Taken for a "+str(nx)+"x"+str(ny)+" Simulation in "+str(size)+" Processors is "+str((end-start)/nt)+" seconds")

 


if __name__ == "__main__":
    args = sys.argv
    nx  = int(args[1])
    ny  = int(args[2])
    nt  = params.nt
    dx  = params.dx
    dy  = params.dy
    dt  = params.dt
    sx  = nx/2    ## Change this for arbitrary source
    sy  = ny/2    ## Change this for arbitrary source
    pad = params.pad
    type = params.type
    checkpoint = params.checkpoint
    output_dir = params.output_dir
    ricker_sigma = params.ricker_sigma
    ricker_shift = params.ricker_shift

    
    size = solve_2d_mpi(nx, ny, nt, sx, sy, dx, dy, dt, pad, type, checkpoint, output_dir)
    

