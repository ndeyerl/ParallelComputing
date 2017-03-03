/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   27 January 2015 */


/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include "vec1d.h"
#include "get_time.h"

/* prototype for Gram-Schmidt routine */
int GramSchmidt1d(vec1d** X, int numvectors);


/* Example routine to test the vec1d "class" */
int main() {
  int i, j, ier;

  /* create some vecs of length 5, and set some entries */
  vec1d *a = vec1dNew(5);
  vec1d *b = vec1dNew(5);
  vec1d *c = vec1dNew(5);
  for (i=0; i<5; i++)  b->data[i] = (i+1)*0.1;
  for (i=0; i<5; i++)  c->data[i] = (i+1);

  /* output to screen */
  printf("writing array of zeros:\n");
  vec1dWrite(a);
  printf("writing array of 0.1, 0.2, 0.3, 0.4, 0.5:\n");
  vec1dWrite(b);
  printf("writing array of 1, 2, 3, 4, 5:\n");
  vec1dWrite(c);

  /* verify that b has length 5 */
  if (b->length != 5)  printf("error: incorrect vector length\n");

  /* update a's data and write to file */
  a->data[0] = 10.0;
  a->data[1] = 15.0;
  a->data[2] = 20.0;
  a->data[3] = 25.0;
  a->data[4] = 30.0;
  vec1dWriteFile(a, "a_data");
  printf("the new file 'a_data' on disk should have entries 10, 15, 20, 25, 30\n\n");
  
  /* access each entry of a and write to screen */
  printf("entries of a, one at a time: should give 10, 15, 20, 25, 30\n");
  for (i=0; i<5; i++) 
    printf("  %g\n", a->data[i]);
  printf("\n");

  /* Test arithmetic operators */
  printf("Testing vector constant, should all be -1\n");
  ier = vec1dConstant(b, -1.0);
  vec1dWrite(b);

  printf("Testing vector copy, should be 1, 2, 3, 4, 5\n");
  ier = vec1dCopy(a, c);
  vec1dWrite(a);

  printf("Testing scalar multiply, should be 5, 10, 15, 20, 25\n");
  ier = vec1dScale(c, 5.0);
  vec1dWrite(c);

  /* create a few vecs of length 10 */
  vec1d* X[5];
  for (i=0; i<5; i++)
    X[i] = vec1dNew(10);
  
  /* fill in the vectors */
  for (i=0; i<10; i++) {
    X[0]->data[i] = 1.0*i;
    X[1]->data[i] = -5.0 + 1.0*i;
    X[2]->data[i] = 2.0 + 2.0*i;
    X[3]->data[i] = 20.0 - 1.0*i;
    X[4]->data[i] = -20.0 + 1.0*i;
  }
  
  /* check the LinearSum routine */
  ier = vec1dLinearSum(X[0], -2.0, X[1], 1.0, X[2]);
  printf("Testing LinearSum, should be all 12:\n");
  vec1dWrite(X[0]);

  /* check the various scalar output routines */
  printf("Testing TwoNorm, should be 2.23607\n");
  printf("  %g\n\n", vec1dTwoNorm(b));
  
  printf("Testing RmsNorm, should be 16.5831\n");
  printf("  %g\n\n", vec1dRmsNorm(c));
  
  printf("Testing MaxNorm, should be 1\n");
  printf("  %g\n\n", vec1dMaxNorm(b));
  
  printf("Testing Min, should be 1\n");
  printf("  %g\n\n", vec1dMin(a));
  
  printf("Testing Max, should be 25\n");
  printf("  %g\n\n", vec1dMax(c));
  
  printf("Testing Dot, should be 275\n");
  printf("  %g\n\n", vec1dDot(a,c));
  
  printf("Testing Linspace, should be 0 1 2 3 4\n");
  vec1d *d = vec1dLinspace(0.0, 4.0, 5);
  vec1dWrite(d);
  
  printf("Testing Random\n");
  vec1d *f = vec1dRandom(5);
  vec1dWrite(f);


  
  /*** performance/validity tests (Gram-Schmidt process) ***/
  int n=1000000;
  printf("Running GramSchmidt1d process\n");
  vec1d** Y = (vec1d**) malloc(5 * sizeof(vec1d*));
  for (i=0; i<5; i++)
    Y[i] = vec1dRandom(n);
  double stime = get_time();
  if (GramSchmidt1d(Y, 5))
    printf("GramSchmidt1d returned error\n");
  double ftime = get_time();
  double rtime = ftime-stime;
  printf("Resulting vectors should be orthonormal:\n");
  for (i=0; i<5; i++)
    for (j=i; j<5; j++)
      printf("  <Y[%i],Y[%i]> = %g\n", i, j, vec1dDot(Y[i],Y[j]));
  printf("\n");
  printf("testing time: %g\n", rtime);


  /* clean up */
  vec1dDestroy(a);
  vec1dDestroy(b);
  vec1dDestroy(c);
  vec1dDestroy(d);
  vec1dDestroy(f);
  for (i=0; i<5; i++)
    vec1dDestroy(X[i]);
  for (i=0; i<5; i++)
    vec1dDestroy(Y[i]);
  free(Y);
  
  return 0;
} /* end main */

