# Ising Model Project, Part 3 (for 440 students)
### MATH 440/540, Spring 2023
##### Code and Final Writeup Due:
11:59 PM Mountain time, May 6, 2023 in GitHub repository
(this is determined by final exam schedule for 11 am MWF classes)

### Collaboration and References
This is intended to be done in groups of 2. A group of 3 will be allowed if someone is unable to find a partner with approval from the instructor. For full credit, all group members are expected to contribute approximately equally (GitHub displays number of lines contributed). One repository is to be submitted per group. If your group discusses this project with another group, you may not copy any code from that collaborator's screen, but may discuss general concepts and share references. You must list every classmate you collaborate with in your references.md file. 

*For this 3rd project part, you are required to have a different partner than you had for the 1st and 2nd projects.* Each portion of the project will require you to find a new partner so that you learn to collaborate with a variety of working styles. For this 3rd project, the 440 and 540 expectations deviate significantly so you and your partner must both be in the 440 section.

### Introduction to Part 3 of the Project

In this part of the project we will build on the part 1 solution to parallelize the process of running a single very large simulation. In this way we will be able to study the behavior of the Ising model under certain conditions with many particles (more like what would happen with a large chunk of a material).   

In particular, you will now work with a simulation that is broken into groups of neighboring rows following a domain decomposition strategy. To progress from one time step to the next, particles on the boundaries of each local subdomain need to have access to the updated states of their neighboring particles on the boundaries of the neighboring processes' local subdomains. Between each time step the ghost regions must be updated. To quantitatively understand the scalability of the code, you will also write strong and weak scalability scripts.

References: 
* The reference most closely related to this part of the project is chapter 2.6 of *"Parallel Algorithms in Computational Science"* by Heermann and Burkitt.
* Some notes on ghost regions by Kjosltad and Snir [https://sites.cs.ucsb.edu/~gilbert/cs140/notes/GhostCells.pdf](https://sites.cs.ucsb.edu/~gilbert/cs140/notes/GhostCells.pdf)


### Deliverables

For the coding portion of this assignment, you will create the following code and produce a markdown writeup with references, all due by 11:59 pm, May 6, 2023 through GitHub. No files other than the files initially provided and the files listed below should be submitted. 
1. **code/Makefile** with rules for basic building and testing examples (provided, you'll make one modification)
2. **code/simulatePar.h** containing struct definition and function prototypes of an individual simulation (provided, won't be modified)
3. **code/simulatePar.c** containing function implementations of an individual parallelized simulation (provided, will be modified)
4. **code/testsProvidedSim.c** provided to test many (but note- not all) parts of your code (provided, you can modify to add extra tests, but don't remove any in the final version submitted)
5. **code/myTest.c** a main program that tests your ghostRegionExchange function (not provided)
6. **code/visualize.py** a Python script provided to visualize results of the provided simulation tests (provided to show results of testsProvidedSim.c)
7. **code/timeSim.c** test case to time your code (provided, won't be modified)
8. **code/weak_scaling.sh** weak scaling test to study scalability (not provided)
9. **code/strong_scaling.sh** strong scaling test to study scalability (provided)
10. **code/visScaling.py** Python script to visualize the results of your strong and weak scaling (provided but you need to modify with your timings copied in)
11. **writeup/writeup.md** markdown file where you will write up your results
12. **writeup/images/_.png** the images referenced in your writeup
13. **writeup/references.md** markdown file where you will write a list of your collaborators and any references you used. 
14. **writeup/slurm-weak.out** text file (you'll rename and move) with results from weak_scaling.sh
15. **writeup/slurm-strong.out** text file (you'll rename and move) with results from strong_scaling.sh

Note: only the files listed above should be pushed to your GitHub repository. No other executable files or output data should be committed.

### Instructions

1. Complete the initializeGridRowBlocks() function to initialize your distributed grid. Check that it can pass the first few tests in testsProvidedSim.c (just comment out the tests for parts you haven't completed yet).

2. Complete the grabUnpaddedData() function to grab the unpadded part of each processes' paddedStateLoc and copy it entry-by-entry into another pre-allocated array. There is a related test in testsProvidedSim that you should use to verify your implementation. 

3. Update the neighborSum() function to calculate the sum of the neighbors of any (row,column) particle. The row and column are a global row and column, and only the rank that actually contains this global row and column in its local (unpadded part of its) subarray will return a correct answer. All other ranks should return a 5 as an indicator of a junk value (as 5 is not a valid value of the sum of neighbors). There are related tests in tests in testsProvidedSum that you should use to verify your implementation, although they are not exhaustive (i.e. you will need to implement additional tests to ensure that your implementation is fully correct). Note: it is okay to review the serial solution from part 1 and part 2 of the project for this, although your code will require significant modifications. 

4. Update the ghostRegionExchange() function to pass the appropriate nPadRows rows from the edge of each process's unpadded subarray to its neighbors and to get the appropriate nPadRows rows from each process's neighboring processes to store in the padded parts of its paddedStateLoc. 

5. Design 3 unit tests for your ghostRegionExchange(), and be sure that they don't just copy the code in your ghostRegionExchange() function (otherwise too likely to let errors pass). Add the implementation of your unit tests into a main program in myTest.c to ensure that your ghost region exchange is working correctly. It must include a clear pass or fail message. Add a rule to the Makefile called myTests to build and run this test with more than one process.

6. Update the hamiltonian() function to calculate the entire Hamiltonian of the whole distributed grid system. Review the formula in project part 1 to remind yourself.  You can assume that the ghost region exchange will be called prior to calling hamiltonian(). First you will fill in the part of the code that aims to calculate the local contribution of each process to the Hamiltonian. Then you will fill in the part of the code that adds these up and calculates the global Hamiltonian of the entire system. Every process should have a copy of the global Hamiltonian upon finishing this function. There is a relevant test in testsProvidedSim. Be careful to ensure that points along the boundaries are not double counted.

7. Update the hamilChange() function to calculate the change in the Hamiltonian due to flipping the sign of a particle's state at any (row,column) particle. The row and column are a global row and column, and only the rank that actually contains this global row and column in its local (unpadded part of its) subarray will return a correct answer. All other ranks should return a 1. There is a relevant test in testsProvidedSim.

8. Complete the timeStep() function with the appropriate indices of each particle within the local padded subarray. A relevant integration test is provided in testsProvidedSim.

9. Complete the simulate() function to calculate the total proportion of particles at each time step throughout the entire global array (will require collective communication) and share this with all processes. All processes should also have access to the global Hamiltonian at each step. 

10. Using the rules in the Makefile and a small job submission script of your own design (does not need to be submitted), run the provided integration tests in testsProvidedSim to ensure that many time steps work in the simulate() function. Run the script in visualize.py as described in the Makefile to visualize some snapshots of the global grid of states, the proportion positive over time, and the Hamiltonian over time. Note this rule is specialized to expect the grid was run with 2 processes, but you could modify this in your own submission script. In the resulting figures, you should expect to see behavior similar to the serial simulations in previous projects.


11. Use the buildSim and runSimExample rules in the Makefile to explore and understand timeSim.c example that times your simulations. Then, run the strong\_scaling.sh script for a strong scaling test on a 2000x2000 grid with 10 repeated simulation tests (for more accurate timing), and the grid is divided over 2, then 4, 8, 16, and 32 processes. 

12. Write a weak scaling study in weak\_scaling.sh that uses the sim executable with P = 2, 4, 8, 16, then 32 processes to simulate a 1000P x 1000 grid, and repeats that simulation 10 times (for more accurate timing results).


13. Modify visScaling.py to include your timings from the strong and weak scaling tests, and run with the appropriate rule in the makefile.


14. Answer the writeup questions below in **writeup/writeup.md** and include all references and collaborators in **writeup/references.md**.

15. **Submit:** You should have been committing code changes regularly. Be sure those changes are pushed to GitHub via the status and log git commands(can even check the website if you aren't sure). LAST  STEP: Be sure to open up your submitted repository on the github.com website and make sure that all the files you expect are there.  I cannot give you credit for things that aren't in your repository. Double check that no binary executable files are in there.


### Writeup Questions (answers must be complete sentences)

1. (6% total) Describe (a, 2%) the initial distribution of the data, (b, 2%) which data are getting sent between processes between time steps, and (c, 2%) which data are used in grabbing the sum of neighbors. You must use specific matrix indices (starting with row 0 and column 0) to describe the data managed and moved between processes.

2. (5%) Describe how you designed your unit test for the ghostRegionExchange() function. Show your work for how you expect the test to behave, and what the expected results are.

3. (6% total) These questions address the weak scalability of your code.

	4a. (2%) Include your weak scaling test's output slurm file in writeup/slurm-weak.out.  
 
 	4b. (2%) Include a plot of your weak scaling timings. What do you observe about your timings?

 	4c. (2%) Fill in this table. What do you observe about your efficiency? What would you expect ideally?

	| n_procs   |      weak test time      |  weak test efficiency |
	|----------|-------------|------|
	| 2 | | |
	| 4 | | |
	| 8 | | |
	| 16 | | |
	| 32 | | |


4. (6% total) These questions address the strong scalability of your code.

	4a. (2%) Include your strong scaling test's output slurm file in writeup/slurm-strong.out.  
 
 	4b. (2%) Include a plot of your strong scaling timings. What do you observe about your timings?

 	4c. (2%) Fill in this table. What do you observe about your efficiency? What would you expect ideally?

	| n_procs   |      strong test time      |  strong test efficiency |
	|----------|-------------|------|
	| 2 | | |
	| 4 | | |
	| 8 | | |
	| 16 | | |
	| 32 | | |


5. (no specific percentage points in the writeup, but this can help you get some points back on code) If you believe your code still has any bugs, you must describe them here so that we can attempt to evaluate your work for partial credit. Which of your test cases shows these bugs (note: you can add an extra test file and makefile rule if needed)? How would you expect the outputs of these test cases change were your code to be fixed? Which lines of code do you believe contain the bug? How have you attempted to fix this portion of code?

6. (2%) references.md must have all collaborators and references listed.


### Grading Rubric
* (10%) Appropriate use of GitHub
* (25%) Code compiles - Note: You will recieve 0% of the 25% if the runTestsProvidedSim rule in the Makefile does not compile on Mio. 
* (20%) Code correctness - will be judged based on provided tests as well as a variety of additional test cases.
* (10%) Scaling scripts are correct and run as expected
* (10%) Coding clarity, organization and comments
* (25%) Writeup and references 
