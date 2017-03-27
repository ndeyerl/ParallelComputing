/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   13 January 2013 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"


// Example routine using the basic 6 MPI functions
int main(int argc, char* argv[]) {

  // local variables
  int ierr, numprocs, myid, n, is, ie, i;
  double *a, *b, stime, ftime, runtime, sum, mysum=0.0;

  // intialize MPI, get total number of processes, this proc's ID
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // set problem parameters
  n = 10000000;

  // root node outputs parallelism information to screen
  if (myid == 0) 
    printf(" starting MPI with %i processes\n",numprocs);

  // determine this proc's interval in overall problem domain
  // (assumes homogeneous system, giving all procs equal work)
  is = floor(1.0*n/numprocs)*myid + 1;
  ie = floor(1.0*n/numprocs)*(myid+1);
  if (myid == numprocs-1)  ie = n;

  // initialize the vectors (only those parts that reside on this proc)
  a = new double[ie-is+1];
  b = new double[ie-is+1];
  for (i=0; i<(ie-is+1); i++) {
    a[i] = 0.001 * (i+is) / n;
    b[i] = 0.001 * (n-i-is) / n;
  }

  // start timer
  stime = MPI_Wtime();

  // perform local portion of dot-product
  for (i=0; i<(ie-is+1); i++)  mysum += a[i]*b[i];

  // root node collects result
  ierr = MPI_Reduce(&mysum, &sum, 1, MPI_DOUBLE, 
		    MPI_SUM, 0, MPI_COMM_WORLD);

  // stop timer
  ftime = MPI_Wtime();
  runtime = ftime - stime;

  // root node outputs compute value and runtime
  if (myid == 0) {
    printf(" dot-product = %.12e\n",sum);
    printf("     runtime = %.12e\n",runtime);
  }

  // free vectors
  delete[] a;
  delete[] b;

  // finalize MPI
  ierr = MPI_Finalize();

} // end main
