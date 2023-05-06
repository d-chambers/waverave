#!/usr/bin/python3
import os

dx = 1
dy = 1
dt = 0.0001
nt = 1000
sx = 0 
sy = 0 
pad = 4
checkpoint = 50
##Sources
s_type = 0

## Ricker
ricker_sigma = 60
ricker_shift = 0.0001 

##Layer
z = [0, 135, 165, 300]
vel = 5000

####OUTPUT
file_path = os.path.abspath(os.path.dirname(__file__))
output_dir = file_path+'/../outputs'
