/* 
 * Nicole Deyerl
 * MATH 6370 Spring 2017
 * File exhibiting errors at compile time by Blaise Blarney and Dan
 * Reynolds, fixed using omp_bug1fix.cpp as a guide. */

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define N           50
#define CHUNKSIZE   5

int main (int argc, char *argv[]) {

  // local variables
  int i, chunk, tid;
  float a[N], b[N], c[N];
  char first_time;
  // Some initializations
  for (i=0; i < N; i++)  a[i] = b[i] = i * 1.0;
  chunk = CHUNKSIZE;

# pragma omp parallel for    \
  shared(a,b,c,chunk)        \
  private(i,tid)             \
  schedule(static,chunk)     \
  firstprivate(first_time)
  for (i=0; i<N; i++) {
    if (first_time == 'y') {
      tid = omp_get_thread_num();
      first_time = 'n';
    }
    c[i] = a[i] + b[i];
    printf("tid = %i,  i = %i,  c[i] = %g\n", tid, i, c[i]);
  }

  return 0;
} // end main

