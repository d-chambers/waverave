# WaveRave!

Term project for Math 540 (CSM, spring 2023)

Derrick Chambers, Peiyao Li, and Ayon Ghosh

<img src="images/waverave.png" width="400" height="350" />



# Introduction

Finite difference methods are a class of numerical methods for solving partial differential equations (PDEs). They are used in a wide variety of applications, including fluid dynamics, heat transfer, and electromagnetism. In this project, we used the finite difference method for solving the wave equation, which is a PDE that describes the propagation of waves. We developed three implementations in different programming languages (C, python, and Julia) all of which used MPI to parallelize the simulations. The following details the domain, implementation tests, and scaling tests.

The acoustic wave equation is given as:

$$
\partial_t^2 p = c^2 \Delta p + s
$$

Where $p$ is the pressure, $c$ is the velocity, $\Delta$ is the Laplacian operator, and $s$ is the source. 

# Simulation Parameters

We selected a simple 2d domain acoustic consisting of two high velocity zones, each 135m thick, surrounding one low velocity zone, 30m thick, and a single explosive source occurring in the center of the domain whose source time function is a 20Hz Ricker wavelet. 

![](images/domain.png)

Grid cells were 1m x 1m, and a temporal sampling was selected such that the CFL criterion (Marfurt, 1984) did not exceed 0.5 for any part of the model. 

$$
V_{max} = \frac{dx}{dt}C_{max}
$$

$$
dt = \frac{dx}{V_{max}}C_{max} = \frac{1}{5000}0.5 = \frac{1}{10000} 
$$


Where $V_{max}$ is the maximum model velocity, $dt$ and $dx$ are the time step and grid discretization, and $C_{max}$ is the maximum CFL criterion.

For simplicity we simply use fixed (reflecting) boundary conditions.




1. (5%) Describe the serial algorithm that you are parallelizing. What problem is it aiming to address? Describe the algorithm with specific sentences or pseudocode. How many operations does the serial algorithm for a particular problem size?

2. (10%) Using detailed sentences or pseudocode, how are you breaking up this algorithm using a distributed memory (MPI) model? Why is this algorithm not embarrassingly parallel? Describe the initial distribution of the data/model, and what data need to be transferred between processes throughout. If it is helpful to introduce new notation, please do so, but be sure to define each index/variable specifically defined at its first use. You may include diagrams as needed. 

3. (5%) Describe your test suite and use its results to make a case that your code is behaving as expected. How do you know your code is right? Are there any limitations to your test suite? 

4. (6%) Describe the aspect of your code you chose to visualize to understand (and quickly check) whether the code is behaving as expected. What did you visualize? How did you do this? How do you know it should look like this? Include some of the figures you produced. 

5. (7%) Describe how you designed your strong and weak scalability tests (dimensions and number of processes). Fill out two tables with the headers (number of procs, weak test time, weak test efficiency) and (number of procs, strong test time, strong test efficiency). Comment on whether your efficiency changes, and make a case for what the biggest contributing factors may be to any slowdowns. 

6. (2%) references.md must have all collaborators and references listed.




# Julia

## Overview and setup

The Julia version of the code is found in the julia folder. To use it in a new account on Mio, the `setup_environment.sh` script needs to be run. This will take a few minutes as it creates a new conda environment with julia installed, then uses julia to instantiate the WaveRave package. Here are some other important files:


`run_jwaverave.jl` - Julia script to run the simulation. Accepts a variety of command line arguments for controlling the simulation. 

`plot_wavefield.jl` - Plots the wavefield (creates a series of pngs at snapshots.)

`submit_weak.sh` - Slurm submission script to get the timing for weak scaling test. Also shows how to run jwaverave with MPI.

`submit_strong.sh` - Slurm submission script to get the timing for strong scaling test



## Testing
There are approximately 30 tests which can be run using the `run_tests.sh` script. These include some simple unit tests as well as a few end2end tests. One test ensures the serial version of the code outputs the same results and the MPI version.

The following gif shows the output of the default serial simulation which uses a homogeneous velocity model.

![](images/julia/animation.gif)

The simulation appears healthy; there is no obvious numerical dispersion and the expected polarity reversals associated with the boundary reflection occur. Since one of the tests compares the output of the MPI version with the serial version we are confident both are accurate.


## Timing



