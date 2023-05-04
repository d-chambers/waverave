#!/usr/bin/python3


import numpy as np
from mpi4py import MPI
import params
import os


 

   
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
          t=np.arange(0, self.dt*self.nt, self.dt)
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
      self.new[self.pad+isx, isy] += f*dt*dt
      

      

   def set_boundaries(self):
       if self.rank == 0:
         self.old[self.pad:2*self.pad, :] = 0
       if self.rank == self.size-1:
         self.old[-2*self.pad:-self.pad,:] = 0

   def ghost_region(self):
      recv_request = []
      send_request = []
      if self.rank > 0:
         receive_previous = self.new[0:self.pad,:]
         recv_request.append(self.comm.Irecv(receive_previous, source = self.rank - 1, tag=1))
         send_previous = self.new[self.pad: 2*self.pad, :]
         send_request.append(self.comm.Isend(send_previous, dest = self.rank - 1, tag=0))

      if self.rank < self.size - 1:
         receive_next = self.new[-self.pad:, :]
         recv_request.append(self.comm.Irecv(receive_next, source = self.rank + 1, tag=0))
         send_next = self.new[-2*self.pad:-self.pad, :]
         send_request.append(self.comm.Isend(send_next, dest = self.rank + 1, tag=1))
      #MPI.Request.waitall(recv_request, recv_status)
      #MPI.Request.waitall(send_request, send_status)
      return recv_request, send_request
   
   def update_solution(self):
      if self.pad == 4:
         
         self.old[4:-4,4:-4] = 2*self.new[4:-4,4:-4]-self.old[4:-4,4:-4]+self.dtdx2*self.v[:,4:-4]**2*( 
                        -1/560 *self.new[0:-8,4:-4] 
                        +8/315 *self.new[1:-7,4:-4]  
                        -1/5   *self.new[2:-6,4:-4]  
                        +8/5   *self.new[3:-5,4:-4]  
                        -205/72*self.new[4:-4,4:-4]  
                        +8/5   *self.new[5:-3,4:-4]  
                        -1/5   *self.new[6:-2,4:-4]  
                        +8/315 *self.new[7:-1,4:-4]  
                        -1/560 *self.new[8:  ,4:-4])+self.dtdy2*self.v[:,4:-4]**2*( 
                        -1/560 *self.new[4:-4,0:-8] 
                        +8/315 *self.new[4:-4,1:-7]  
                        -1/5   *self.new[4:-4,2:-6]  
                        +8/5   *self.new[4:-4,3:-5]  
                        -205/72*self.new[4:-4,4:-4]  
                        +8/5   *self.new[4:-4,5:-3]  
                        -1/5   *self.new[4:-4,6:-2]  
                        +8/315 *self.new[4:-4,7:-1]  
                        -1/560 *self.new[4:-4,8:  ])
      else:
         raise RuntimeError("Method not yet implemented for padding size other than 4")
      
   
   def write(self, time):
       np.save(os.path.dirname(os.path.dirname(__file__))+'/outputs/wavefield_'+str(time).zfill(5)+'_'+str(self.rank).zfill(5), self.new[self.pad:-self.pad,:])
                    