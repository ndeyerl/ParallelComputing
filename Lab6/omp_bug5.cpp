/* 
 * Nicole Deyerl
 * MATH 6370 Spring 2017
 * Failed compilation bug by Blaise Blarney and Dan Reynolds, fixed 
 * by changing how program works to model omp_orphan.cpp by DR. */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define VECLEN 100

// global variables
float a[VECLEN], b[VECLEN], sum;

// external function
void dotprod() {

  // local variables 
  int i, tid;

  tid = omp_get_thread_num();
# pragma omp for reduction(+:sum)
  for (i=0; i < VECLEN; i++) {
    sum = sum + (a[i]*b[i]);
    printf("  tid = %i,  i = %i\n", tid, i);
  }
}


// main routine
int main (int argc, char *argv[]) {

  // local variables
  int i;

  // initialize values
  for (i=0; i<VECLEN; i++)   a[i] = b[i] = 1.0 * i;
  sum = 0.0;

  // begin parallel region, with orphaned call to dotprod()
# pragma omp parallel
  dotprod();

  // output result
  printf("Sum = %g\n", sum);

  return 0;
} // end main
