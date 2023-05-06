#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from mpi4py import MPI
import params 
import simulation




def plot_2d(array):
    fig = plt.imshow(array, cmap="seismic")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    
    return fig




def plot_vel_one_rank(nx, ny, nt, dx, dt, output_dir,n):

    #output_dir = os.path.dirname(os.path.dirname(__file__))+'outputs'

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    global_wave = np.zeros((nx,ny))
    start = 0
    for i in range(n):
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


def plot_wave_one_rank(nx, ny, nt, dx, dt, checkpoint, output_dir,n):
    if checkpoint != 0:
       time = np.arange(checkpoint,nt,checkpoint)
       if nt not in time:
          time = np.append(time,nt)
    else:
        time =[nt]

    #output_dir = os.path.dirname(os.path.dirname(__file__))+'outputs'




    for t in time:
        global_wave = np.zeros((nx,ny))
        start = 0
        for i in range(n):
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
   n = int(args[3])
   nt = params.nt 
   dt = params.dt
   dx = params.dx
   dy = params.dy
   checkpoint = params.checkpoint
   output_dir = params.output_dir
   os.system("mkdir "+output_dir+"/plots")
   plot_wave_one_rank(nx, ny, nt, dx, dt, checkpoint, output_dir,n)
   plot_vel_one_rank(nx, ny, nt, dx, dt, output_dir,n)

