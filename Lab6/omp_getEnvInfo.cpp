/* FILE: omp_getEnvInfo.cpp
   Nicole Deyerl
   Math6370 Spring 2017
   OpenMP file written to get info about the OpenMP environment, modified   
   from Blaise barney's Hello World OpenMP example which was updated to 
   c++ by Dan Reynolds.*/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {

  // local variables
  int nprocs, nthreads, mthreads, ifpar, ifdpar, ifnpar;

  // Fork a team of threads giving them their own copies of variables
# pragma omp parallel private(nthreads, tid)
  {

    // Only master thread does this
    if (tid == 0) {
      nprocs = omp_get_num_procs();
      printf("Number of processors = %i\n", nprocs);

      nthreads = omp_get_num_threads();
      printf("Number of threads = %i\n", nthreads);

      mthreads = omp_get_max_threads();
      printf("Max number of threads available = %i\n", mthreads);

      ifpar = omp_in_parallel();
      printf("In parallel region = %i\n",ifpar);

      ifdpar = omp_get_dynamic();
      printf("Dynamic threads enabled = %i\n",ifdpar);

      ifnpar = omp_get_nested;
      printf("Nested parallelism supported = %i\n",idnpar);
    }

  } // All threads join master thread and disband
  
  return 0;
}  // end main
