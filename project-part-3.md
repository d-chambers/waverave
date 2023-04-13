## Final Project - Distributed Memory Parallelism
### MATH 540, Spring 2023
##### Code and Final Writeup Due:
11:59 PM Mountain time, May 6, 2023 in GitHub repository
(this is determined by final exam schedule for 11 am MWF classes)

### Collaboration and References
This is intended to be done in groups of 2. A group of 3 will be allowed if someone is unable to find a partner with approval from the instructor. For full credit, all group members are expected to contribute approximately equally (GitHub displays number of lines contributed). One repository is to be submitted per group. If your group discusses this project with another group, you may not copy any code from that collaborator's screen, but may discuss general concepts and share references. You must list every classmate you collaborate with in your references.md file. 

*For this 3rd project part, you are required to have a different partner than you had for the 1st and 2nd projects.* Each portion of the project will require you to find a new partner so that you learn to collaborate with a variety of working styles. For this 3rd project, the 440 and 540 expectations deviate significantly so you and your partner must both be in the 540 section.


### Introduction to Part 3 of the Project

In this part of the project you will practice parallel algorithm design, implementation of a distributed memory parallel algorithm, and testing to ensure that the algorithm behaves as expected (correctness and scalability). 


### Deliverables

For the coding portion of this assignment, you will create the following code and produce a markdown writeup with references, all due by 11:59 pm, May 6, 2023 through GitHub. 
0. **project-part-3.md** this project description
1. **code/README.md** a readme markdown file where you will list every file you have in your code subdirectory and its general description (1-2 lines per file)
2. All .c and .h files used to implement your algorithm
3. All main programs used for scalability testing and for running unit/integration tests.
4. **code/test.sh** script to run all of your unit and integration tests
5. **code/weak_scaling.sh** weak scaling test to study scalability
6. **code/strong_scaling.sh** strong scaling test to study scalability 
7. **code/scriptsForVisualization** Any scripts you wrote to produce your visualizations.
8. **writeup/writeup.md** markdown file where you will write up your results
9. **writeup/images/_.png** the images referenced in your writeup
10. **writeup/references.md** markdown file where you will write a list of your collaborators and any references you used. 
11. **writeup/slurm-weak.out** text file (you'll rename and move) with results from weak_scaling.sh
12. **writeup/slurm-strong.out** text file (you'll rename and move) with results from strong_scaling.sh

Note: only the files (or file types) listed above should be pushed to your GitHub repository. No executable files or output data should be committed.

### Instructions

1. Select an algorithm that you would like to parallelize using MPI. This algorithm cannot be embarrassingly parallel to qualify as an appropriate choice for this project. Please come to office hours or email the instructor to discuss your algorithm if you are unsure whether it would qualify for this project. 

2. Answer the writeup questions below in **writeup/writeup.md** relevant to the parallel algorithm design, and include all references and collaborators in **writeup/references.md**.

3. Decide on the interfaces you'd like to use for structs and functions in your code. Write a suite of unit tests and integration tests (at least 8 unit tests and 4 integration tests) to build confidence that your code is working as expected. One of your tests must compare serial and parallel results to show that they are equal up to expected numerical precision error.

4. Implement your parallel algorithm. Be sure to break your code into smaller functions and use structures as appropriate (and this will also help with your testing suite).

5. Create a job submission script for running your test suite in parallel on a small number of cores (4-12 cores). Use your test suite to help guide your debugging process. Once your code is debugged, answer all writeup questions below in **writeup/writeup.md** relevant to the parallel algorithm testing. 

6. Create a script (in Python, Matlab, VisIt, VTK or any other visualization program you prefer) to visualize some kind of result from your method that you can use to help explain what is going on in your code. 
For example, if you implement a parallel matrix LU decomposition, you might show images of the matrix values before factorization, images of the L and U factor values of each factor after factorization, and a plot of the factor products compared to the original matrix (to show that they are nearly identical if the code is correct).
Answer all writeup questions below in **writeup/writeup.md** to show this visualization.

7. Design a strong scalability test and a weak scalability test, write job submission scripts for these, and run them. Answer all writeup questions below in **writeup/writeup.md** about your scalability. 

8. **Submit:** You should have been committing code changes regularly. Be sure those changes are pushed to GitHub via the status and log git commands(can even check the website if you aren't sure). LAST  STEP: Be sure to open up your submitted repository on the github.com website and make sure that all the files you expect are there.  I cannot give you credit for things that aren't in your repository. Double check that no binary executable files are in there.


### Writeup Questions (answers must be complete sentences)


1. (5%) Describe the serial algorithm that you are parallelizing. What problem is it aiming to address? Describe the algorithm with specific sentences or pseudocode. How many operations does the serial algorithm for a particular problem size?

2. (10%) Using detailed sentences or pseudocode, how are you breaking up this algorithm using a distributed memory (MPI) model? Why is this algorithm not embarrassingly parallel? Describe the initial distribution of the data/model, and what data need to be transferred between processes throughout. If it is helpful to introduce new notation, please do so, but be sure to define each index/variable specifically defined at its first use. You may include diagrams as needed. 

3. (5%) Describe your test suite and use its results to make a case that your code is behaving as expected. How do you know your code is right? Are there any limitations to your test suite? 

4. (6%) Describe the aspect of your code you chose to visualize to understand (and quickly check) whether the code is be behaving as expected. What did you visualize? How did you do this? How do you know it should look like this? Include some of the figures you produced. 

5. (7%) Describe how you designed your strong and weak scalability tests (dimensions and number of processes). Fill out two tables with the headers (number of procs, weak test time, weak test efficiency) and (number of procs, strong test time, strong test efficiency). Comment on whether your efficiency changes, and make a case for what the biggest contributing factors may be to any slowdowns. 

6. (2%) references.md must have all collaborators and references listed.


### Grading Rubric
* (5%) Appropriate use of GitHub
* (30%) Code compiles, runs, and reproduces correct results on Mio with fully reproducible scripts.
* (20%) Testing suite that you design (complete and follows best practices)
* (10%) Coding clarity, organization and comments
* (35%) Writeup and references 
