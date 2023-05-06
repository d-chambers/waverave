from mpi4py import MPI
import numpy as np
import sys
import os
from src import simulation, solver
from src import params

def test_source(nx, ny, dx, dy, dt,sx, sy, pad):

    comm=MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    
    par2d = simulation.simul_par(nx, ny , dx, dy, comm, rank, size)
    par2d.size_allocate()
    rank_source, isx, isy = par2d.source(sx,sy,pad)
    
    if rank == rank_source:
       if rank_source != 2:
           print("Test 1 Failed: Source is not injected at correct Rank which is 2 in the test case")
       else:
           print("Test 1 Passed: Source injected at current Rank")
    comm.Barrier()


def test_decomposition(nx, ny, dx, dy, dt):


    comm=MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    par2d = simulation.simul_par(nx, ny, dx, dy, comm, rank, size)

    par2d.size_allocate()

    if par2d.lnx == 20:
        print("Test 2 of Decomposition Size Passed: For Evenly divisible domain and MPI size, every process gets same number of chunks")
    else:
        print("Test 2 of Decomposition Size Failed: Domain Decomposition does not work")
    comm.Barrier()


def test_velocity(nx, ny, dx, dy, dt, pad, output_dir):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    par2d = simulation.simul_par(nx, ny , dx, dy, comm, rank, size)
    par2d.size_allocate()
    v = np.zeros((par2d.lnx, par2d.lny))
    v[:,:] = rank
    #print(np.shape(v), rank)

    np.save(output_dir+'/vel_'+str(rank).zfill(5),v)
    
    file = os.path.abspath(os.path.dirname(__file__))+'/outputs/vel_'+str(rank).zfill(5)+'.npy'
    wavefield = solver.wavefield(par2d.lnx, par2d.lny, pad, dx, dy, dt, comm, rank, size, par2d.method)
    wavefield.vel(file)
    wavefield.new[:,:] = rank

    if wavefield.method == 0:
        diff = wavefield.new[pad:-pad,:]- wavefield.v[:,:]
        if np.any(diff):
            print("Test 3 Failed: Wrong velocity model being read from "+str(rank))
        else:
            print("Test 3 Passed: Correct velocity model being read from rank "+str(rank))
    elif wavefield.method == 1:
        diff = wavefield.new[:,pad:-pad]- wavefield.v[:,:]
        if np.any(diff):
            print("Test 3 Failed: Wrong velocity model being read from "+str(rank))
        else:
            print("Test 3 Passed: Correct velocity model being read from rank "+str(rank))
    comm.Barrier()


def test_region(nx, ny, dx, dy,dt, pad):
    comm=MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    par2d = simulation.simul_par(nx, ny, dx, dy, comm, rank, size)

    par2d.size_allocate()

    wavefield = solver.wavefield(par2d.lnx, par2d.lny, pad, dx, dy, dt, comm, rank, size, par2d.method)

    wavefield.new[:,:] = rank

    wavefield.set_exchange()

    recv_request, send_request = wavefield.ghost_region()
    recv_status=[]
    send_status=[]
    MPI.Request.waitall(recv_request, recv_status)
    MPI.Request.waitall(send_request, send_status)
    comm.Barrier()
    
    test_score = 0
    if nx>=ny:
        cmp_array = np.zeros((pad,ny))
        if rank > 0:
           cmp_array[:,:] = rank-1
           cmp_array=wavefield.new[0:pad,:] - cmp_array
           if np.any(cmp_array):  
              test_score +=1
        if rank < size-1:
            cmp_array[:,:] = rank+1
            cmp_array = wavefield.new[-pad:,:] - cmp_array
            if np.any(cmp_array):   
                test_score +=1
    else:
        cmp_array = np.zeros((nx,pad))
        if rank > 0:
           cmp_array[:,:] = rank-1
           cmp_array=wavefield.new[:,0:pad] - cmp_array
           if np.any(cmp_array):  
              test_score +=1
        if rank < size-1:
            cmp_array[:,:] = rank+1
            cmp_array = wavefield.new[:,-pad:] - cmp_array
            if np.any(cmp_array):   
                test_score +=1
    
    
    if test_score != 0:
            print("Test 4 Failed: Ghost Region Exchange Not Working in Rank "+str(rank))
    else:
            print("Test 4 Passed: Ghost Region Exchange Working in Rank "+str(rank))
    


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
    output_dir = os.path.abspath(os.path.dirname(__file__))+"/outputs"
    
    test_source(nx, ny, dx, dy, dt, sx, sy, pad)
    print("")
    print("")
    test_decomposition(nx, ny, dx, dy, dt)
    print("")
    print("")
    test_velocity(nx, ny, dx, dy, dt, pad, output_dir)
    print("")
    print("")
    test_region(nx, ny, dx, dy,dt, pad)
