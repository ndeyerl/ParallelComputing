/* 
 * Nicole Deyerl
 * MATH 6370 Spring 2017
 * Buggy file by Blaise Blarney and Dan Reynolds, fixed using Dan 
 * Reynolds' hints.  */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {

  // local variables
  int nthreads, i, tid;
  float total;

  // Spawn parallel region 
# pragma omp parallel private(tid) //made tid private
  {
    // Obtain thread number
    tid = omp_get_thread_num();
    
    // Only master thread does this
    if (tid == 0) {
      nthreads = omp_get_num_threads();
      printf("Number of threads = %i\n", nthreads);
    }

    // Everyone does this
    printf("Thread %i is starting...\n", tid);

#   pragma omp barrier

    // do some work
    total = 0.0;
#   pragma omp for schedule(dynamic,10) private(total) //made total private
    for (i=0; i<1000000; i++)  total = total + i*1.0;

    printf("Thread %i is done! Total = %g\n", tid, total);

  } // End of parallel region

  return 0;
}  // end main

