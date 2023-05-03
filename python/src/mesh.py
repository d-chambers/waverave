#!/usr/bin/python3

import sys
import domain 

def main(argv,dx,dy,vel):
    if len(argv) != 3:
        print("Only give nx ny as input")
        return None
    else:
         nx = float(sys.argv[1])
         nz = float(sys.argv[2])
         xmin = 0
         ymin = 0
         xmax = nx*dx
         ymax = nz*dy
         extent = [xmin, ymin, xmax, ymax]
         
         domain.make_homogeneous(extent, dx, dy, vel)


if __name__ ==  '__main__':
    print(sys.argv)
    main(sys.argv, 1, 1, 4000)
