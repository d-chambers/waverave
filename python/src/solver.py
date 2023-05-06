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
   def __init__ (self, local_nx, local_ny, pad, dx, dy, dt, comm, rank, size, file, method):
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
      self.method = method
      v = np.load(file)
      lnx = np.shape(v)[0]
      lny = np.shape(v)[1]
      if lnx == self.lnx and lny == self.lny:
          self.v = v
      else:
           raise IndexError("Rerun Mesher: Velocity Grid of rank "+str(self.rank)+ " has different dimensions from expected")
      if self.method == 0:
         self.old = np.zeros((self.lnx+2*self.pad, self.lny))
         self.new = np.zeros((self.lnx+2*self.pad, self.lny)) 
      else:
         self.old = np.zeros((self.lny+2*self.pad, self.lnx))
         self.new = np.zeros(( self.lny+2*self.pad, self.lnx,))  
         self.old=self.old.T
         self.new=self.new.T        
   
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
         self.old[self.b0[0]:self.b0[1],self.b0[2]:self.b0[3]] = 0
       if self.rank == self.size-1:
         self.old[self.b1[0]:self.b1[1],self.b1[2]:self.b1[3]]  = 0

   def set_exchange(self):
       if self.method == 0:
          r_previous = [0, self.pad, 0, self.lny] 
          s_previous = [self.pad, 2*self.pad, 0, self.lny]
          r_next = [self.lnx + self.pad, self.lnx +2*self.pad, 0, self.lny]
          s_next = [self.lnx, self.lnx+self.pad, 0, self.lny]
          boundary_0 = [self.pad, 2*self.pad, 0, self.lny]
          boundary_1 = [self.lnx, self.lnx+self.pad,0, self.lny]
          v =[0,self.lnx,self.pad,self.lny-self.pad]
       else:
           r_previous = [0, self.lnx, 0, self.pad] 
           s_previous = [0, self.lnx, self.pad, 2*self.pad]
           r_next = [0, self.lnx, self.lny + self.pad, self.lny+2*self.pad]
           s_next = [0, self.lnx, self.lny, self.lny+self.pad]
           boundary_0 = [0, self.lnx, self.pad, 2*self.pad ]
           boundary_1 = [0, self.lnx, self.lny, self.lny+self.pad]
           v = [self.pad, self.lnx-self.pad,  0, self.lny]
       self.rp = r_previous
       self.sp = s_previous
       self.rn = r_next
       self.sn = s_next
       self.b0 = boundary_0
       self.b1 = boundary_1
       self.ve = v
       print(self.rp,self.rn,self.sp,self.sn)

   def ghost_region(self):
      recv_request=[]
      send_request=[] 
      #receive_previous = np.zeros((self.rp[1]-self.rp[0], self.rp[3]-self.rp[2]))
      #receive_next = np.zeros((self.rn[1]-self.rn[0], self.rn[3]-self.rn[2]))
      if self.rank > 0:
         receive_previous = self.new[self.rp[0]:self.rp[1],self.rp[2]: self.rp[3]]
         recv_request.append(self.comm.Irecv(receive_previous, source = self.rank - 1, tag=1))
         #self.new[self.rp[0]:self.rp[1],self.rp[2]: self.rp[3]] = receive_previous
         send_previous = self.new[self.sp[0]:self.sp[1],self.sp[2]: self.sp[3]]
         send_request.append(self.comm.Isend(send_previous, dest = self.rank - 1, tag=0))

      if self.rank < self.size - 1:
         
         receive_next = self.new[self.rn[0]:self.rn[1],self.rn[2]: self.rn[3]]
         recv_request.append(self.comm.Irecv(receive_next, source = self.rank + 1, tag=0))
         #self.new[self.rn[0]:self.rn[1],self.rn[2]: self.rn[3]] = receive_next
         send_next = self.new[self.sn[0]:self.sn[1],self.sn[2]: self.sn[3]]
         print(np.max(send_next))
         send_request.append(self.comm.Isend(send_next, dest = self.rank + 1, tag=1))
      #MPI.Request.waitall(recv_request, recv_status)
      #MPI.Request.waitall(send_request, send_status)
      return recv_request, send_request
   
   def update_solution(self):
      if self.pad == 4:
         self.old[4:-4,4:-4] = 2*self.new[4:-4,4:-4]-self.old[4:-4,4:-4]+self.dtdx2*self.v[self.ve[0]:self.ve[1],self.ve[2]:self.ve[3]]**2*( 
                        -1/560 *self.new[0:-8,4:-4] 
                        +8/315 *self.new[1:-7,4:-4]  
                        -1/5   *self.new[2:-6,4:-4]  
                        +8/5   *self.new[3:-5,4:-4]  
                        -205/72*self.new[4:-4,4:-4]  
                        +8/5   *self.new[5:-3,4:-4]  
                        -1/5   *self.new[6:-2,4:-4]  
                        +8/315 *self.new[7:-1,4:-4]  
                        -1/560 *self.new[8:  ,4:-4])+self.dtdy2*self.v[self.ve[0]:self.ve[1],self.ve[2]:self.ve[3]]**2*( 
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
      
   
   def write(self, time, output_dir):
       if self.method==0:
          np.save(output_dir+'/wavefield_'+str(time).zfill(5)+'_'+str(self.rank).zfill(5), self.new[self.pad:-self.pad,:])
       else:
          np.save(output_dir+'/wavefield_'+str(time).zfill(5)+'_'+str(self.rank).zfill(5), self.new[:,self.pad:-self.pad])
                    