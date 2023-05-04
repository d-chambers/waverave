#!/usr/bin/python3
class simul_par(object):
   def __init__(self, nx, ny, dx, dy, comm, rank, size):
       """
       nx: Dimension in x direction
       ny: Dimension in y direction
       nt: Number of time steps in simulation
       dx: Grid Spacing in x direction
       dy: Grid Spacing in y direction
       dt: Time step
       sx: Source Location in x direction
       sy: Source Location in y direction
       F:  Forcing array at each time step
       comm: MPI_Communicator
       """

       self.nx = nx
       self.ny = ny
       self.dx = dx
       self.dy = dy
       self.comm = comm
       self.rank = rank
       self.size = size
       


   def size_allocate(self):
      """
      Gives Memory Sizes for a given rank and total size
      nx: Size of the total data to be parallelized(Here number of rows)
      size: Number of Processes in MPi
      rank: Rank of the processor it's working
      """
      chunk_size = int(self.nx / self.size)
      if self.rank < self.nx % self.size:
         chunk_size+=1
         start= self.rank*chunk_size
      else:
         start = self.rank*chunk_size + self.nx % self.size
      
      end = start+chunk_size

      self.start = start
      self.chunk_size = int(chunk_size)
      self.end = end

   def source(self, sx, sy):
      """ 
      Function to inject source in a proper location in a MPI case
      rank_source: rank of the domain where source should be injected
      isx: x location of the source in the local wavefield
      isy: y location of the source in the local wavefield 
      """

      isx = int(sx / self.dx)
      isy = int(sy / self.dy)
      if (isx-self.start) < self.chunk_size and (isx-self.start)>=0:
        rank_source = self.rank
        isx = isx-self.start
        print("isx="+str(isx))
        print(self.chunk_size)
        isy = isy
        return rank_source, isx, isy
      else: 
           return self.size+1, isx, isy 