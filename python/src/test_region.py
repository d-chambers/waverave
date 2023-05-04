import solver
import params
import simulation
import os
from mpi4py import MPI

#args = sys.argv
nx  = 15
ny  = 15
nt  = params.nt
dx  = params.dx
dy  = params.dy
dt  = params.dt
sx  = nx/2    ## Change this for arbitrary source
sy  = ny/2    ## Change this for arbitrary source
pad = params.pad
type = params.type
checkpoint = params.checkpoint

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

par2d = simulation.simul_par(nx, ny , dx, dy, comm, rank, size)
par2d.size_allocate()

file = os.path.dirname(os.path.dirname(__file__))+'/outputs/vel_'+str(rank).zfill(5)+'.npy'

wavefield = solver.wavefield(par2d.chunk_size, par2d.ny, pad, dx, dy, dt, comm, rank, size, file)

wavefield.new[:,:] = wavefield.rank
wavefield.old[:,:] = wavefield.rank



wavefield.write(0)

print(wavefield.new, rank)
rs=[]
ss=[]
r,s = wavefield.ghost_region()
MPI.Request.waitall(r,rs)
MPI.Request.waitall(s,ss)
print(wavefield.new, rank, "ghosted")

wavefield.write(100)

wavefield.set_boundaries()
#print(wavefield.new, rank)
print(wavefield.old,rank, "biundaried")