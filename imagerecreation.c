#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "src/pgmio.h"
#include "src/log.h"
#include "src/inputfuncs.h"

#define MAXITER  1500
#define PRINTFREQ  200

#define ndims 2
#define TRUE 1
#define FALSE 0

const char *FILENAME = "edge_files/edgenew192x128.pgm";

float boundaryval(int i, int m);
float get_avg_pixel_val(int M, int N, float arr[M][N], int not_include_halos);

int main (int argc, char **argv)
{

	/* ----------------------------------------------------- */
	/* MPI Initiation */
	MPI_Init (NULL, NULL);
	int rank, size, left, right, up, down;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
	MPI_Request request;
	MPI_Status status;
	int TAG = 0;

	/* ----------------------------------------------------- */
	/* Define and set Input Parameters from the Command line */
	int M, N, Mp, Np;
	int i, j, iter, sawtooth_iter;
	float val;
	int num_iters= MAXITER, print_freq= PRINTFREQ;
	double stop_delta = 0.001;
	char *filename = strdup(FILENAME);
	
	if (argc >1) {
		set_input_params(argc, argv, &num_iters, &print_freq, filename, &stop_delta);
	} 
	print_input_params(rank, num_iters, print_freq, filename, stop_delta);	

	/* ----------------------------------------------------- */
	/* Timer Variables */
	double starttime, endtime;
	double parallel_starttime, parallel_endtime;
	double file_read_start, file_read_end, file_write_start, file_write_end;
	starttime = MPI_Wtime();

	/* ----------------------------------------------------- */	
	/* Read the input file Size */
	pgmsize(filename, &M,&N);

	/* ----------------------------------------------------- */
	/* Stopping Criterion variables
	 * stop_loop flag will be used to stop the iterations when criteria is fulfilled
	 * check_freq is the frequency of running this check.
	 */
	float current_delta, max_delta = 0.0;
	float max_delta_all_procs=0.0;
	int check_freq = num_iters/5; 
	int stop_loop = 0;

	/* ----------------------------------------------------- */
	/* Cartesian topology
	   2d topology, cyclic on the first dimension, and not on the second one
	*/
	int comm2d,direction_1d,direction_2d,disp;
	int dims[ndims], period[ndims], coords[ndims];
	int coords_x[size], coords_y[size];
	int reorder;
	int x0,x1,y0,y1;
	dims[0] = 0;
	dims[1] = 0;
	coords[0] = 0;
	coords[1] = 0;
	period[0] = TRUE;    // TRUE, Cyclic horizontal
	period[1] = FALSE;   // not cyclic vertical
	reorder = FALSE;
	direction_1d = 0;       // shift along the first index
	direction_2d = 1;	//shift along the 2nd index
	disp = 1;            // Shift by 1

        MPI_Dims_create(size,ndims,dims);
        MPI_Cart_create(comm,ndims,dims,period,reorder,&comm2d);
        MPI_Comm_rank(comm2d,&rank);
        MPI_Cart_shift(comm2d,direction_1d,disp,&left,&right);
	MPI_Cart_shift(comm2d,direction_2d,disp,&up,&down);
	MPI_Cart_coords(comm2d , rank, ndims, coords);

	log_it(rank, "The dimensions of my Cartesian Coordinates are: X=%d, Y=%d\n",dims[0],dims[1]);
	debug("Cartesian coordinates: My rank=%d -- left=%d, right=%d, up=%d, down=%d, mpi_null=%d\n",
			rank, left, right, up, down, MPI_PROC_NULL); 
	debug("Cartesian grid breakdown: my_rank=%d, my_X=%d, my_Y=%d\n",
			rank, coords[0], coords[1]);

	/* ----------------------------------------------------- */
	/* Setting the array size values based on the filename and number of procs 
	 * Does not handle cases when M or N is not divisible by the dimensions. 
	/* Using static array allocation
	 */
	Mp = M/dims[0];
	Np = N/dims[1];
	float masterbuf[M][N], buf[Mp][Np];
	float old[Mp+2][Np+2], new[Mp+2][Np+2], edge[Mp+2][Np+2];

	/* ----------------------------------------------------- */
	/* MPI derived data type
	 * Vector of size MpxNp for master to buffer comm
	 * Vector of size Mpx1 for up down halo swaps. Note stride length includes Halos
	 */
	MPI_Datatype vectorMpxNp, updownVector;
	MPI_Type_vector(Mp, Np, Np*dims[1], MPI_FLOAT, &vectorMpxNp);
	MPI_Type_commit( &vectorMpxNp );
	MPI_Type_vector(Mp, 1, Np+2, MPI_FLOAT, &updownVector);
	MPI_Type_commit( &updownVector);

	/* Avg pixel value calcs */
	float avg_pixel, avg_pixel_sum=0.0;


	/* ----------------------------------------------------- */
	/* File IO being done only on rank zero.
	 * This file is then distributed across to the other procs by Send/ Recvs
	 */
	if (0==rank){
		log_it(rank, "Reading <%s>\n", filename);
		file_read_start = MPI_Wtime();
  		pgmread(filename, masterbuf, M, N);	
		file_read_end = MPI_Wtime();
	}


	/* ----------------------------------------------------- */
	/* So far all the items have been serial in nature, including file IO 
	 * From here on, we start measuring the timing of the parallel code
	 */
	parallel_starttime = MPI_Wtime();

	/* All procs send their cartesian coordinates to rank 0 to keep track */
	MPI_Gather(&(coords[0]), 1, MPI_INT, &coords_x[0], 1, MPI_INT, 0, comm2d);
	MPI_Gather(&(coords[1]), 1, MPI_INT, &coords_y[0], 1, MPI_INT, 0, comm2d);
	for (i=0;i<size && rank==0;i++){
		debug("coords_x=%d, coords_y=%d from rank=%d\n", coords_x[i], coords_y[i], i);
	}

	if (rank==0){
		/* Distributing the masterbuf accross to smaller buffers in every proc
		 * For rank 0, we can simply copy the initial portion of array into buf array of rank 0
		 * For the other ranks, we do a synchronous send of vectorMpxNp
		 * The starting index of masterbuf is based on 2d decomp coordinates
		 */
		for (i=0;i<Mp;i++) {
			for(j=0;j<Np;j++) {
				buf[i][j] = masterbuf[i][j];
			}
		}
		for (i=1;i<size;i++) {
			MPI_Ssend(&(masterbuf[coords_x[i]*Mp][coords_y[i]*Np]), 1, vectorMpxNp, i, TAG, comm2d);
		}
	}
	else {
		MPI_Recv(&buf[0][0], Mp*Np, MPI_FLOAT, 0,TAG, comm2d, &status);
	}

	/* ----------------------------------------------------- */
	/* Set up the initial calculation arrays
	 * old array starts out completely white (255.0) for all values
	 * It includes fixed boundary conditions for top and bottom edges
	 * Edge is a copy of buf array
	 * combining the three loops into a single loop for efficiency
	 */
	for (i=0;i<Mp+2;i++)
	{
		for (j=0;j<Np+2;j++)
		{
			old[i][j]=255.0;
			if ((i>=1 && i<=Mp) && (j>=1 && j<=Np)) {
                                edge[i][j]=buf[i-1][j-1];
                        }
		}
		/* Set fixed boundary conditions on the top and bottom edges
		 * compute sawtooth value. 
		 * Assumption 1: Sawtooth formula should apply to the overall image, not individual process
		 * 		This means the iterator i passed to boundaryval func should be 
		 *		the value of i in the overall picture
		 * Assumption 2: This only needs to apply to top and bottom edges of the actual picture
				and not every part of the image. So we check for ranks with up and down 
		 * 		direction being null 
		 */
		if (i>=1 && i<=Mp) {
			sawtooth_iter = i + coords[0]*Mp;
			val = boundaryval(sawtooth_iter, M);
			if (up==MPI_PROC_NULL) {
				old[i][0]   = 255.0*val;
			}
			if (down==MPI_PROC_NULL) {
				old[i][Np+1] = 255.0*(1.0-val);
			}
		}

	}


	for (iter=1;iter<=num_iters && stop_loop ==0; iter++)
	{
		if(iter%print_freq==0)
		{
			log_it(rank, "Iteration %d, rank = %d\n", iter, rank);

			/* Avg Pixel value calculation 
			 * Calc the avg on each proc and then reduce it on rank 0.
			*/
			avg_pixel = get_avg_pixel_val(Mp+2, Np+2, old, TRUE);
			debug ("Rank = %d, Avg pixel value=%f\n", rank, avg_pixel);
			MPI_Reduce(&avg_pixel, &avg_pixel_sum, 1, MPI_FLOAT, MPI_SUM, 0, comm2d);
			log_it(rank, "Avg pixel value = %f\n", avg_pixel_sum/size);
			
		}

		/* ----------------------------------------------------- */
		/* Halo Swaps Implementation - Asynchronous Send & Recvs */

		/* Implement periodic boundary conditions on left and right sides 
		 * Left/Right swaps are easy in C as we are working with 
		 * contiguous memory set. 
		 */
		
		MPI_Issend(&old[Mp][1], Np, MPI_FLOAT, right, TAG, comm2d, &request);
		MPI_Recv(&old[0][1], Np, MPI_FLOAT, left, TAG, comm2d, &status);
		MPI_Wait(&request, &status);

		MPI_Issend(&old[1][1], Np, MPI_FLOAT, left, TAG, comm2d, &request);
		MPI_Recv(&old[Mp+1][1], Np, MPI_FLOAT, right, TAG, comm2d, &status);
		MPI_Wait(&request, &status);
		
		/* Implement periodic boundary conditions on up and down sides 
		 * Up/Down swaps are trickier in C, as we are working with non
		 * contiguous memory set. Need to use a vector which will have the right stride length.
		 */
		
		MPI_Issend(&old[1][Np], 1, updownVector, down, TAG, comm2d, &request);
		MPI_Recv(&old[1][0], 1, updownVector, up, TAG, comm2d, &status);
		MPI_Wait(&request, &status);
	
		MPI_Issend(&old[1][1], 1, updownVector, up, TAG, comm2d, &request);
		MPI_Recv(&old[1][Np+1], 1, updownVector, down, TAG, comm2d, &status);
		MPI_Wait(&request, &status);
		
		for (i=1;i<Mp+1;i++)
		{
			for (j=1;j<Np+1;j++)
			{
				new[i][j]=0.25*(old[i-1][j]
						+old[i+1][j]
						+old[i][j-1]
						+old[i][j+1]
			      			- edge[i][j]);
				current_delta = fabs(new[i][j] - old[i][j]);
				if (current_delta > max_delta) max_delta = current_delta;
			}
		}
		/* Reset the array for the next iteration */
	      	
		for (i=1;i<Mp+1;i++)
		{
		  for (j=1;j<Np+1;j++)
		    {
		      old[i][j]=new[i][j];
		    }
		}

		/* ----------------------------------------------------- */
		/* Stopping criteria check. 
		 * If criteria is met criteria, stop_loop flag will be updated on every proc and
		 * all of them will exit the loop.
		 */
		if (iter!=0 && iter%check_freq == 0) {
			MPI_Allreduce(&max_delta, &max_delta_all_procs, 1, MPI_FLOAT, MPI_MAX, comm2d);
			if (max_delta_all_procs <= stop_delta) {
				log_it(rank, "Iter = %d, max_delta=%f is lower than stop_delta=%f. \n"
						"Stopping further calcs.\n", 	
						iter, max_delta_all_procs, stop_delta);
				stop_loop=1;
			}
			max_delta_all_procs = 0.0;
			max_delta = 0.0;
		}
		
	}//end n iterations or after stop_loop flag is set.

	
	log_it(rank, "Finished %d iterations\n", iter-1);


	for (i=1;i<Mp+1;i++)
	{
		for (j=1;j<Np+1;j++)
		{
			buf[i-1][j-1]=old[i][j];
		}
	}

	/* ----------------------------------------------------- */
	/* This is the gather implementation for this process
	 * Once the image processing is complete, each process
	 * except rank 0, sends its buf array to rank 0, which then
	 * combines the image back into masterbuf array using the 
	 * same indexing as was used to distribute the image first
	 */
	if (rank != 0) {
		MPI_Ssend(&(buf[0][0]), Mp*Np, MPI_FLOAT, 0, TAG, comm2d);
	}
	else{
		for (i=0;i<Mp;i++) {
			for(j=0;j<Np;j++) {
				masterbuf[i][j] = buf[i][j];
			}
		}
		for (i=1;i<size;i++) {
			MPI_Recv(&(masterbuf[(coords_x[i])*Mp][(coords_y[i])*Np]), 1, vectorMpxNp, i, TAG, comm2d, &status);
		}
	}

	/* End of the parallel section */
	parallel_endtime = MPI_Wtime();
		
	if (rank ==0){
		sprintf(filename,"imagenew%dx%d.pgm",M, N);
		log_it(rank, "Writting masterbuf to output file\n");
		file_write_start = MPI_Wtime();
		pgmwrite(filename, masterbuf, M, N);
		file_write_end = MPI_Wtime();
	}

	endtime = MPI_Wtime();
	log_it(rank, "Run time split (seconds)\n"
			"\tTotal Run time= %f\n"
			"\tParallel Run time= %f\n"
			"\tFile IO time = %f\n"
			,endtime-starttime, parallel_endtime - parallel_starttime
			,(file_read_end-file_read_start)+(file_write_end-file_write_start)
			);
	
	MPI_Finalize();

	return 0;
} 

float boundaryval(int i, int m)
{
  float val;

  val = 2.0*((float)(i-1))/((float)(m-1));
  if (i >= m/2+1) val = 2.0-val;
  
  return val;
}

/* Calcs the avg pixel value for a given array.
 * Simply sums the vals of the array and divides by the total count to get the avg
 * If not_include_halos is true , we will remove halo effects from this calc 
 *  -- start and end points changes by 1 iter
 */
float get_avg_pixel_val(int M, int N, float arr[M][N], int not_include_halos)
{
	float res, pixel_sum =0;
	int i, j;
	int halo_effect;

	if (not_include_halos == TRUE) { 
		halo_effect = 1;
	} else {
		halo_effect = 0; 
	}

	for (i= (0+halo_effect) ; i < (M-halo_effect) ; i++){
		for(j= (0+halo_effect) ; j < (N-halo_effect) ; j++) {
			pixel_sum += arr[i][j];
		}
	}
	res = pixel_sum/(M*N);
	return res;
}
