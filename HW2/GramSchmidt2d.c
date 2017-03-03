/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   10 February 2017 */


/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include "vec2d.h"


/* tolerance for linear independence checks */
#define TOL 1e-12

/* Gram-Schmidt process */
int GramSchmidt2d(vec2d** X, int numvectors) {
  int i, j;
  double Xnorm;
  
  /* check that there is work to do */
  if (numvectors < 1)  return 0;
  
  /* normalize first vector */
  Xnorm = vec2dTwoNorm(X[0]);
  if (fabs(Xnorm) < TOL) {
    fprintf(stderr,"GramSchmidt2d error: vectors are linearly-dependent!\n");
    return 1;
  } else {
    vec2dScale(X[0], 1.0/Xnorm);
  }
  
  /* iterate over remaining vectors, performing Gram-Schmidt process */
  for (i=1; i<numvectors; i++) {

    /* subtract off portions in directions of existing basis vectors */
    for (j=0; j<i; j++)
      vec2dLinearSum(X[i], 1.0, X[i], -vec2dDot(X[i],X[j]), X[j]);
    
    /* normalize vector, checking for linear dependence */
    Xnorm = vec2dTwoNorm(X[i]);
    if (fabs(Xnorm) < TOL) {
      fprintf(stderr,"GramSchmidt2d error: vectors are linearly-dependent!\n");
      return 1;
    } else {
      vec2dScale(X[i], 1.0/Xnorm);
    }
  }

  /* return success */
  return 0;
}
