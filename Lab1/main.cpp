/* Description: 

   Reads input file 'input.txt' to obtain vector length n.
   Creates two vectors of length n, u and v, and performs the 
   dot-product of these two vectors, outputing the result both 
   to screen and to the file 'output.txt'

   Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   11 January 2013
 */

// Inclusions
#include <stdlib.h>   // new, delete
#include <time.h>     // clock()
#include <stdio.h>    // printf()


int main(int argc, char* argv[]) {

  // declarations
  int i, n;
  double *u, *v, sum;
  FILE *FID;

  /* read problem parameters from input file (should be in this order):
        nx - number of nodes in x-direction
	ny - number of nodes in y-direction
	nt - number of time steps
	tstop - final time (will stop at nt or stop, whichever is 1st)
	c - wave speed
	dtoutput - time frequency for outputting solutions */
  FID=fopen("input.txt","r");
  fscanf(FID,"n = %i,\n", &n);
  fclose(FID);

  printf("  vector length = %i\n", n);

  // allocate arrays u and v
  u = new double[n];
  v = new double[n];

  // fill in vectors u and v
  for (i=0; i<n; i++) {
    u[i] = i+1;
    v[i] = n-i;
  }

  // compute dot-product
  sum = 0.0;
  for (i=0; i<n; i++)   sum += u[i]*v[i];

  // output computed value and runtime
  printf("  dot-product = %g\n",sum);

  // output length and dot-product to file 'output.txt'
  FID=fopen("output.txt","w");
  fprintf(FID, " n = %i\n", n);
  fprintf(FID, " u dot v = %g\n", sum);
  fclose(FID);

  // delete vectors u and v
  delete[] u;
  delete[] v;

  return 0;
} // end main
