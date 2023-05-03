#!/usr/bin/python3


import numpy as np
from mpi4py import MPI
import params
import os



# class simul_par():
#    def __init__(self, nx, ny, dx, dy, comm, rank, size):
#        """
#        nx: Dimension in x direction
#        ny: Dimension in y direction
#        nt: Number of time steps in simulation
#        dx: Grid Spacing in x direction
#        dy: Grid Spacing in y direction
#        dt: Time step
#        sx: Source Location in x direction
#        sy: Source Location in y direction
#        F:  Forcing array at each time step
#        comm: MPI_Communicator
#        """

#        self.nx = nx
#        self.ny = ny
#        self.dx = dx
#        self.dy = dy
#        self.comm = comm
#        self.rank = rank
#        self.size = size


#    def size_allocate(self):
#       """
#       Gives Memory Sizes for a given rank and total size
#       nx: Size of the total data to be parallelized(Here number of rows)
#       size: Number of Processes in MPi
#       rank: Rank of the processor it's working
#       """
#       chunk_size = int(self.nx / self.size)
#       if self.rank < self.nx % self.size:
#          chunk_size+=1
#          start= self.rank*chunk_size
#       else:
#          start = self.rank*chunk_size + self.nx % self.size
      
#       end = start+chunk_size

#       self.start = start
#       self.chunk_size = chunk_size
#       self.end = end

#    def source(self, sx, sy):
#       """ 
#       Function to inject source in a proper location in a MPI case
#       rank_source: rank of the domain where source should be injected
#       isx: x location of the source in the local wavefield
#       isy: y location of the source in the local wavefield 
#       """

#       isx = int(sx / self.dx)
#       isy = int(sy / self.dy)

#       if (isx-self.start) < self.chunk_size and (isx-self.start)>=0:
#          rank_source = self.rank
#          isx = isx-self.start
#          isy = isy
#          return rank_source, isx, isy


 

   
class Ricker(object):
   def __init__(self, nt, dt):
      """
      Creates Ricker Wavelet Class Instance
      ss = sigma for Ricker wavelet
      t0 = time shift of the wavelet
      nt = total number of time step in the simulation
      dt = time step 
      """
      self.ss = params.ricker_sigma
      self.t0 = params.ricker_shift
      self.dt = dt
      self.nt = nt
      
   def source_time(self):
          F=np.zeros(self.nt)
          t=np.linspace(0, self.dt*self.nt, self.dt)
          F = (1-((t*self.dt - self.t0)/self.ss)**2)*np.exp(-(t*self.dt - self.t0)**2/(2*self.ss**2))
          return F
      

class wavefield(object):
   def __init__ (self, local_nx, local_ny, pad, dx, dy, dt, comm, rank, size, file):
      """
      lnx = dimension in x direction in local array
      lny = dimension in y direction in local 
      pad = padding dimensions in ghost regions
      stencil = stencil array showing coefficients for various displacement
      rank = rank of the processor in which it is working
      size = size of the MPI pool
      comm = MPI Communicator
      """
      self.lnx = local_nx
      self.lny = local_ny
      self.pad   = pad
      self.rank = rank
      self.size = size
      self.comm = comm
      self.dt = dt
      self.dx = dx
      self.dy = dy
      self.dtdx2 = (dt/float(dx))**2
      self.dtdy2 = (dt/float(dy))**2
      v = np.load(file)
      lnx = np.shape(v)[0]
      lny = np.shape(v)[1]
      if lnx == self.lnx and lny == self.lny:
          self.v = v
      else:
           raise IndexError("Rerun Mesher: Velocity Grid of rank "+str(self.rank)+ " has different dimensions from expected")
      self.old = np.zeros((self.lnx+2*self.pad, self.lny))
      self.new = np.zeros((self.lnx+2*self.pad, self.lny)) 
   
   # def vel(self, file):
   #    v = np.load(file)
   #    lnx = np.shape(v)[0]
   #    lny = np.shape(v)[1]
   #    if lnx == self.lnx and lny == self.lny:
   #        self.v = np.zeros((self.nx, self.ny))
   #    else:
   #         raise IndexError("Rerun Mesher: Velocity Grid of rank "+str(self.rank)+ " has different dimensions from expected")
      
   def inject_source(self, isx, isy, f, dt):
      self.new[isx, isy] += f*dt*dt

      

   def set_boundaries(self):
       if self.rank == 0:
         self.old[self.pad:2*self.pad, :] = 0
       if self.rank == 1:
         self.old[-2*self.pad:self.pad,:] = 0
       self.old[:,:self.pad] = 0
       self.old[:,-self.pad:] = 0

   def ghost_region(self):
      recv_request = []
      send_request = []
      if self.rank > 0:
         bottom_ghost_next = self.new[0:self.pad,:]
         recv_request.append(self.comm.Irecv(bottom_ghost_next, source = self.rank - 1, tag=1))
         #self.new[0:4, :] = bottom_ghost_next
         top_ghost_next = self.new[self.pad: 2*self.pad, :]
         send_request.append(self.comm.Isend(top_ghost_next, dest = self.rank - 1, tag=0))

      if self.rank < self.size - 1:
         top_ghost_previous = self.new[-self.pad:, :]
         recv_request.append(self.comm.Irecv(top_ghost_previous, source = self.rank + 1, tag=0))

         bottom_ghost_previous = self.new[-8:-4, :].copy()

         send_request.append(self.comm.Isend(bottom_ghost_previous, dest = self.rank + 1, tag=1))
      return recv_request, send_request
   
   def update_solution(self):
      if self.pad == 4:
         self.old[4:4+self.lnx,4:self.lny-4] = 2*self.new[4:4+self.lnx,4:self.lny-4]-self.old[4:4+self.lnx,4:self.lny-4]+self.dtdx2*self.v[:,4:self.lny-4]**2*( 
                        -1/560 *self.new[0:4+self.lnx-4,4:self.lny-4] 
                        +8/315 *self.new[1:4+self.lnx-3,4:self.lny-4]  
                        -1/5   *self.new[2:4+self.lnx-2,4:self.lny-4]  
                        +8/5   *self.new[3:4+self.lnx-1,4:self.lny-4]  
                        -205/72*self.new[4:4+self.lnx,4:self.lny-4]  
                        +8/5   *self.new[5:4+self.lnx+1,4:self.lny-4]  
                        -1/5   *self.new[6:4+self.lnx+2,4:self.lny-4]  
                        +8/315 *self.new[7:4+self.lnx+3,4:self.lny-4]  
                        -1/560 *self.new[8:4+self.lnx+4  ,4:self.lny-4])+self.dtdy2*self.v[:,4:self.lny-4]**2*( 
                        -1/560 *self.new[4:4+self.lnx,0:self.lny-8] 
                        +8/315 *self.new[4:4+self.lnx,1:self.lny-7]  
                        -1/5   *self.new[4:4+self.lnx,2:self.lny-6]  
                        +8/5   *self.new[4:4+self.lnx,3:self.lny-5]  
                        -205/72*self.new[4:4+self.lnx,4:self.lny-4]  
                        +8/5   *self.new[4:4+self.lnx,5:self.lny-3]  
                        -1/5   *self.new[4:4+self.lnx,6:self.lny-2]  
                        +8/315 *self.new[4:4+self.lnx,7:self.lny-1]  
                        -1/560 *self.new[4:4+self.lnx,8:self.lny  ])
      else:
         raise RuntimeError("Method not yet implemented for padding size other than 4")
      
   
   def write(self, time):
       np.save(os.path.dirname(os.path.dirname(__file__))+'/outputs/wavefield_'+str(time).zfill(5)+'_'+str(self.rank).zfill(5), self.new)
                    