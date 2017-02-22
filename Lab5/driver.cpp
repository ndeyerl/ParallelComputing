/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   10 February 2015 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <omp.h>

// Prototypes
void vector_sum(int, double, double*, double, double*, double*);
void vector_scale(int, double, double*, double*);
void vector_product(int, double*, double*, double*);
void vector_pow(int, double*, double, double*);
double vector_rmsnorm(int, double*);


/* Allocates and fills the vectors x, y, z.
   Calls a set of vector routines.
   Writes the total time taken for all vector routines. */
int main(int argc, char* argv[]) {

  // local variables
  int n, i;
  double *x, *y, *z, a, b;
  double stime, ftime;

  // set problem parameters
  n = 10000000;

  // allocate data arrays
  x = new double[n];
  y = new double[n];
  z = new double[n];

  // fill in vectors x and y, scalars a and b
  for (i=0; i<n; i++) {
    x[i] = 1.0*(i+1)/n;
    y[i] = 1.0*(n-i-1)/n;
  }
  a = -2.0;
  b = 0.2;

  // start timer -- Note: we use the OpenMP timer omp_get_wtime()
  // instead of clock(), since omp_get_wtime() is independent of the
  // number of threads, but clock() adds the time on each thread together
  stime = omp_get_wtime();

  // call the vector routines a number of times
  for (i=0; i<20; i++) {
    vector_sum(n, b, x, a, y, z);
    vector_scale(n, b, z, y);
    vector_product(n, y, z, x);
    vector_pow(n, x, b, y);
    a = vector_rmsnorm(n, y);
  }

  // stop timer
  ftime = omp_get_wtime();

  // output result and solution time
  std::cout << " Run using " << omp_get_max_threads() << " threads\n";
  std::cout << " Final rms norm = " << std::setprecision(16) << a << std::endl;
  std::cout << " Total run time = " << ftime - stime << std::endl;

  // free data arrays
  delete[] x;
  delete[] y;
  delete[] z;

  return 0;
} // end main
