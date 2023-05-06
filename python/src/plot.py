#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from mpi4py import MPI
import params 
import simulation

def plot_wave(nx, ny, nt, dx, dt, checkpoint):
    if checkpoint != 0:
       time = np.arange(checkpoint,nt,checkpoint)
       if nt not in time:
          time = np.append(time,nt)
    else:
        time =[nt]

    output_dir = os.path.dirname(os.path.dirname(__file__))+'outputs'

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    for t in time:
        vel = simulation.simul_par(nx, ny, dx, dy, comm, rank, size)
        vel.size_allocate()


        
        
        y_min = vel.start * vel.dx
        y_max = vel.end * vel.dx
        x_min = 0
        x_max = ny * vel.dy
        extent = [x_min, x_max, nx*dx-y_min, nx*dx-y_max, vel.ny, vel.chunk_size]
        if rank != 0:
            comm.send(extent, dest = 0)
        if rank == 0:
            fig, ax = plt.subplots(figsize=(10,10))
            wave = np.load(output_dir+'/wavefield_'+str(t).zfill(5)+'_'+str(0).zfill(5)+'.npy')

            for s in range(1,size):
                x = np.linspace(extent[0], extent[1], extent[4])
                y = np.linspace(extent[2], extent[3], extent[5])
                X,Y = np.meshgrid(x,y)
                ax.pcolormesh(X,Y,wave)
                extent = comm.recv(source=s)
                wave = np.load(output_dir+'/wavefield_'+str(t).zfill(5)+'_'+str(s).zfill(5)+'.npy')
            x = np.linspace(extent[0], extent[1], extent[4])
            y = np.linspace(extent[2], extent[3], extent[5])
            X,Y = np.meshgrid(x,y)
            ax.pcolormesh(X,Y,wave)
            plt.xlabel("X (m)")
            plt.ylabel("Y (m)")
            plt.title(" time = "+str(t*dt)+" seconds")
            fig.savefig(output_dir+"/plots/Wavefield_"+str(t).zfill(5)+"png")


def plot_vel(nx, ny, nt, dx, dt):

    output_dir = os.path.dirname(os.path.dirname(__file__))+'outputs'

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    vel = simulation.simul_par(nx, ny, dx, dy, comm,rank, size)
    vel.size_allocate()

    


    y_min = vel.start * vel.dx
    y_max = vel.end * vel.dx
    x_min = 0
    x_max = ny * vel.dy
    extent = [x_min, x_max, nx*dx-y_min, nx*dx-y_max, vel.ny, vel.chunk_size]
    if rank != 0:
        comm.send(extent, dest = 0)
    if rank == 0:
        fig, ax = plt.subplots()
        wave = np.load(output_dir+'/vel_'+str(0).zfill(5)+'.npy')
        print(wave)
        for s in range(1,size):
            x = np.linspace(extent[0], extent[1], extent[4])
            y = np.linspace(extent[2], extent[3], extent[5])
            X,Y = np.meshgrid(x,y)
            ax.pcolormesh(X,Y,wave)
            extent = comm.recv(source=s)
            wave = np.load(output_dir+'/vel_'+str(s).zfill(5)+'.npy')
            print(wave)
        x = np.linspace(extent[0], extent[1], extent[4])
        y = np.linspace(extent[2], extent[3], extent[5])
        X,Y = np.meshgrid(x,y)
        ax.pcolormesh(X,Y,wave)
        plt.xlabel("X (m)")
        plt.ylabel("Y (m)")
        plt.title(" velocity model")
        fig.savefig(output_dir+"/plots/velocity.png")


def plot_2d(array):
    fig = plt.imshow(array, cmap="seismic")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    
    return fig




def plot_vel_one_rank(nx, ny, nt, dx, dt, output_dir):

    #output_dir = os.path.dirname(os.path.dirname(__file__))+'outputs'

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    global_wave = np.zeros((nx,ny))
    start = 0
    for i in range(size):
        wave = np.load(output_dir+'/vel_'+str(i).zfill(5)+'.npy')
        if np.shape(wave)[1]==ny:
            x = np.shape(wave)[0]
            end = start +x
            global_wave[start:end,:] = wave
            start = end
        else:
            y = np.shape(wave)[1]
            end = start + y
            global_wave[:,start:end] = wave
            start = end


    plot_2d(global_wave)
    plt.savefig(output_dir+"/plots/velocity.png")


def plot_wave_one_rank(nx, ny, nt, dx, dt, checkpoint, output_dir):
    if checkpoint != 0:
       time = np.arange(checkpoint,nt,checkpoint)
       if nt not in time:
          time = np.append(time,nt)
    else:
        time =[nt]

    #output_dir = os.path.dirname(os.path.dirname(__file__))+'outputs'

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    for t in time:
        global_wave = np.zeros((nx,ny))
        start = 0
        for i in range(size):
            wave = np.load(output_dir+'/wavefield_'+str(t).zfill(5)+'_'+str(i).zfill(5)+'.npy')
            if np.shape(wave)[1]==ny and np.shape(wave)[0]!=nx:
               x = np.shape(wave)[0]
               end = start +x
               global_wave[start:end,:] = wave
               start = end
            else:
                y = np.shape(wave)[1]
                end = start + y
                global_wave[:,start:end] = wave
                start = end
            
        plot_2d(global_wave)
        plt.savefig(output_dir+"/plots/wavefield_"+str(t).zfill(5)+"png")





        




if __name__ == "__main__":
   args = sys.argv
   nx  = int(args[1])
   ny  = int(args[2])
   nt = params.nt 
   dt = params.dt
   dx = params.dx
   dy = params.dy
   checkpoint = params.checkpoint
   output_dir = params.output_dir
   plot_wave_one_rank(nx, ny, nt, dx, dt, checkpoint, output_dir)
   plot_vel_one_rank(nx, ny, nt, dx, dt, output_dir)

