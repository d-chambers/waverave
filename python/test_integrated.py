import numpy as np
import sys

def get_diff(nx, ny, n1, n2):
    global_wave_1 = np.zeros((nx,ny))
    start = 0
    t = 250
    for i in range(n1):
        wave = np.load('./outputs_with_'+str(n1)+'_process'+'/wavefield_'+str(t).zfill(5)+'_'+str(i).zfill(5)+'.npy')
        if np.shape(wave)[1]==ny and np.shape(wave)[0]!=nx:
            x = np.shape(wave)[0]
            end = start +x
            global_wave_1[start:end,:] = wave
            start = end
        else:
            y = np.shape(wave)[1]
            end = start + y
            global_wave_1[:,start:end] = wave
            start = end

    global_wave_2 = np.zeros((nx,ny))
    start = 0
    for i in range(n2):
        wave = np.load('./outputs_with_'+str(n2)+'_process'+'/wavefield_'+str(t).zfill(5)+'_'+str(i).zfill(5)+'.npy')
        if np.shape(wave)[1]==ny and np.shape(wave)[0]!=nx:
            x = np.shape(wave)[0]
            end = start +x
            global_wave_2[start:end,:] = wave
            start = end
        else:
            y = np.shape(wave)[1]
            end = start + y
            global_wave_2[:,start:end] = wave
            start = end

    diff = global_wave_1-global_wave_2
    max_diff = np.max(np.abs(diff))

    print("Maximum Difference between the two wavefields at 250th iteration is"+str(max_diff))



if __name__ == "__main__":
    args = sys.argv
    nx  = int(args[1])
    ny  = int(args[2])
    n1 = int(args[3])
    n2 = int(args[4])

    get_diff(nx,ny,n1,n2)
