# parallelImageProc_MPI
Parallel Image processing using MPI
*******************************************************************************
MPP Coursework
Date: December 2015
*******************************************************************************
*******************************************************************************

This folder contains the image recreation program built for MPP Coursework.
To build:
1. Ensure that PrgEnv-pgi environment is setup
2. make

To build with DEBUG mode set ON (to see additional print statements)
1. Edit the Makefile - set debug = <blank> (unset DNDEBUG flag)
2. make

To run:
1. Standalone: ./imagerecreation
2. Using MPI (running on p processors): 
	mpiexec -np <p> ./imagerecreation
3. Passing non default input parameters: 
	mpmiexec -np <p> ./imagerecreation -n 10000 -p 1000 -f edge_files/edgenew256x192.pgm -s 0.01
4. Get help on input parameters:
	./imagerecreation --help

To submit jobs to morar (p Processors):
1. qsub -q "morar1+2" -pe mpi <p> imagerecreation.sge

*******************************************************************************
*******************************************************************************

Folders:
1. edge_files folder contains all the edge files, including the default one which is hard coded in the script
2. src contains all the header files which supply additional functions to the main

