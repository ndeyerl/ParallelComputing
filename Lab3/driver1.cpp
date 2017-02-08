/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   13 January 2013 */


// Inclusions
#include <stdlib.h>   // new, delete
#include <stdio.h>    // printf()
#include <time.h>     // clock()
#include <math.h>     // pow()


// Prototypes
void vector_sum(int, int, int, double***, double***, double***);
void vector_scale(int, int, int, double, double***, double***);
double vector_infnorm(int, int, int, double***);
double vector_dotprod(int, int, int, double***, double***);



/* Allocates and fills the multi-dimensional vectors x, y, z 
   with random numbers.
   Calls a set of vector routines.
   Writes the total time taken for all vector routines. */
int main(int argc, char* argv[]) {

  // local variables
  int l, m, n, i, j, k;
  double ***x, ***y, ***z, a, runtime;
  clock_t stime, ftime;

  // set problem parameters
  l = 100;   // length of first index
  m = 100;   // length of second index
  n = 100;   // length of third index

  // allocate 3D data arrays
  x = new double**[l];
  y = new double**[l];
  z = new double**[l];
  for (i=0; i<l; i++)  x[i] = new double*[m];
  for (i=0; i<l; i++)  y[i] = new double*[m];
  for (i=0; i<l; i++)  z[i] = new double*[m];
  for (i=0; i<l; i++)  
    for (j=0; j<m; j++)  x[i][j] = new double[n];
  for (i=0; i<l; i++)  
    for (j=0; j<m; j++)  y[i][j] = new double[n];
  for (i=0; i<l; i++)  
    for (j=0; j<m; j++)  z[i][j] = new double[n];

  // fill in vectors x and y
  for (i=0; i<l; i++) 
    for (j=0; j<m; j++) 
      for (k=0; k<n; k++) 
	x[i][j][k] = random() / (pow(2.0,31.0) - 1.0);
  for (i=0; i<l; i++)  
    for (j=0; j<m; j++) 
      for (k=0; k<n; k++) 
	y[i][j][k] = random() / (pow(2.0,31.0) - 1.0);

  // start timer
  stime = clock();

  // call the vector routines a number of times
  for (i=0; i<100; i++) {
    vector_sum(l, m, n, x, y, z);
    vector_scale(l, m, n, 2.0, x, z);
    a = vector_infnorm(l, m, n, z);
    a = vector_dotprod(l, m, n, x, y);
  }

  // stop timer
  ftime = clock();
  runtime = ((double) (ftime - stime))/CLOCKS_PER_SEC;

  // output solution time
  printf(" Result from computation = %.12e\n",a);
  printf(" Total run time = %.12e\n",runtime);

  // free data arrays
  for (i=0; i<l; i++) {
    for (j=0; j<m; j++) {
      delete[] x[i][j];
      delete[] y[i][j];
      delete[] z[i][j];
    }
    delete[] x[i];
    delete[] y[i];
    delete[] z[i];
  }
  delete[] x;
  delete[] y;
  delete[] z;

} // end main

