/* Nicole Deyerl
   SMU Mathematics
   Math 4370/6370
   3 March 2017 */


/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include "vec2d_b.h"
#include "get_time.h"

/* prototype for Gram-Schmidt routine */
int GramSchmidt2d_b(vec2d_b** X, int numvectors);


/* Example routine to test the vec2d_b "class" */
int main() {
  int i, j, ier;

  /* create some vecs of length 2x3, and set some entries */
  vec2d_b *a = vec2d_bNew(2,3);
  vec2d_b *b = vec2d_bNew(2,3);
  vec2d_b *c = vec2d_bNew(2,3);
  for (i=0; i<2; i++){
	  for (j=0; j<3; j++){
		  b->data[i*b->length2 + j] = (j+1)*0.1 + i*3.0*0.1;
	  }
  }
  for (i=0; i<2; i++){
	  for (j=0; j<3; j++){
		  c->data[i*c->length2 + j] = (j+1) + i*3.0;
	  }
  }

  /* output to screen */
  printf("writing array of zeros:\n");
  vec2d_bWrite(a);
  printf("writing array of 0.1, 0.2, 0.3, 0.4, 0.5, 0.6:\n");
  vec2d_bWrite(b);
  printf("writing array of 1, 2, 3, 4, 5, 6:\n");
  vec2d_bWrite(c);

  /* verify that b has length 2x3 */
  if ((b->length1!=2) || (b->length2!=3))  printf("error: incorrect vector length\n");
  
  /* update a's data and write to file */
  a->data[0] = 10.0;
  a->data[1] = 15.0;
  a->data[2] = 20.0;
  a->data[3] = 25.0;
  a->data[4] = 30.0;
  a->data[5] = 35.0;
  vec2d_bWriteFile(a, "a_data");
  printf("the new file 'a_data' on disk should have entries 10, 15, 20, 25, 30, 35\n\n");
  
  /* access each entry of a and write to screen */
  printf("entries of a, one at a time: should give 10, 15, 20, 25, 30, 35\n");
  for (i=0; i<2; i++){
	  for (j=0; j<3; j++){
		  printf("  %g\n", a->data[i*a->length2 + j]);
	  }
  }
  printf("\n");

  /* Test arithmetic operators */
  printf("Testing vector constant, should all be -1\n");
  ier = vec2d_bConstant(b, -1.0);
  vec2d_bWrite(b);

  printf("Testing vector copy, should be 1, 2, 3, 4, 5, 6\n");
  ier = vec2d_bCopy(a, c);
  vec2d_bWrite(a);

  printf("Testing scalar multiply, should be 5, 10, 15, 20, 25, 30\n");
  ier = vec2d_bScale(c, 5.0);
  vec2d_bWrite(c);

  /* create a few (5) vecs of length 2x2 */
  vec2d_b* X[5];
  for (i=0; i<5; i++){
	  X[i] = vec2d_bNew(2,2);
  }

  
  /* fill in the vectors */
  for (i=0; i<2; i++){
	  for (j=0; j<2; j++){
		  X[0]->data[i*X[0]->length2 + j] = 1.0*i;
		  X[1]->data[i*X[1]->length2 + j] = 1.0*j;
		  X[2]->data[i*X[2]->length2 + j] = 1.0*i + 1.0*j;
		  X[3]->data[i*X[3]->length2 + j] = 1.0*i + 5.0*j;
		  X[4]->data[i*X[4]->length2 + j] = 5.0*i + 1.0*j;
	  }
  }

  
  /* check the LinearSum routine */
  ier = vec2d_bLinearSum(X[0], -2.0, X[1], 1.0, X[2]);
  printf("Testing LinearSum, should be [0 -1; 1 0]:\n");
  vec2d_bWrite(X[0]);

  /* check the various scalar output routines */
  printf("Testing TwoNorm, should be 2.4495\n");
  printf("  %g\n\n", vec2d_bTwoNorm(b));
  
  printf("Testing RmsNorm, should be 19.4722\n");
  printf("  %g\n\n", vec2d_bRmsNorm(c));
  
  printf("Testing MaxNorm, should be 1\n");
  printf("  %g\n\n", vec2d_bMaxNorm(b));
  
  printf("Testing Min, should be 1\n");
  printf("  %g\n\n", vec2d_bMin(a));
  
  printf("Testing Max, should be 30\n");
  printf("  %g\n\n", vec2d_bMax(c));
  
  printf("Testing Dot, should be 455\n");
  printf("  %g\n\n", vec2d_bDot(a,c));
  
  printf("Testing Linspace, should be 0 1 2 3 4 5 6 7 8 9 10 11\n");
  vec2d_b *d = vec2d_bLinspace(0.0, 11.0, 4, 3);
  vec2d_bWrite(d);
  
  printf("Testing Random\n");
  vec2d_b *f = vec2d_bRandom(5, 2);
  vec2d_bWrite(f);



  /*** performance/validity tests (Gram-Schmidt process) ***/

  printf("Running GramSchmidt2d_b process\n");

  vec2d_b** Y1 = (vec2d_b**) malloc(5 * sizeof(vec2d_b*));
  for (i=0; i<5; i++)
    Y1[i] = vec2d_bRandom(10000, 1000);
  double stime = get_time();
  if (GramSchmidt2d_b(Y1, 5))
    printf("GramSchmidt2d_b returned error\n");
  double ftime = get_time();
  double rtime = ftime-stime;
  printf("Resulting vectors should be orthonormal:\n");
  for (i=0; i<5; i++)
    for (j=i; j<5; j++)
      printf("  <Y1[%i],Y1[%i]> = %g\n", i, j, vec2d_bDot(Y1[i],Y1[j]));
  printf("\n");
  printf("dimensions: (10000, 1000)\n");
  printf("testing time: %g\n", rtime);

  vec2d_b** Y2 = (vec2d_b**) malloc(5 * sizeof(vec2d_b*));
  for (i=0; i<5; i++)
    Y2[i] = vec2d_bRandom(1000, 10000);
  stime = get_time();
  if (GramSchmidt2d_b(Y2, 5))
    printf("GramSchmidt2d_b returned error\n");
  ftime = get_time();
  rtime = ftime-stime;
  printf("Resulting vectors should be orthonormal:\n");
  for (i=0; i<5; i++)
    for (j=i; j<5; j++)
      printf("  <Y1[%i],Y1[%i]> = %g\n", i, j, vec2d_bDot(Y2[i],Y2[j]));
  printf("\n");
  printf("dimensions: (1000, 10000)\n");
  printf("testing time: %g\n", rtime);
  
  vec2d_b** Y3 = (vec2d_b**) malloc(5 * sizeof(vec2d_b*));
  for (i=0; i<5; i++)
    Y3[i] = vec2d_bRandom(100, 100000);
  stime = get_time();
  if (GramSchmidt2d_b(Y3, 5))
    printf("GramSchmidt2d_b returned error\n");
  ftime = get_time();
  rtime = ftime-stime;
  printf("Resulting vectors should be orthonormal:\n");
  for (i=0; i<5; i++)
    for (j=i; j<5; j++)
      printf("  <Y1[%i],Y1[%i]> = %g\n", i, j, vec2d_bDot(Y3[i],Y3[j]));
  printf("\n");
  printf("dimensions: (100, 100000)\n");
  printf("testing time: %g\n", rtime);  
  
  vec2d_b** Y4 = (vec2d_b**) malloc(5 * sizeof(vec2d_b*));
  for (i=0; i<5; i++)
    Y4[i] = vec2d_bRandom(10, 1000000);
  stime = get_time();
  if (GramSchmidt2d_b(Y4, 5))
    printf("GramSchmidt2d_b returned error\n");
  ftime = get_time();
  rtime = ftime-stime;
  printf("Resulting vectors should be orthonormal:\n");
  for (i=0; i<5; i++)
    for (j=i; j<5; j++)
      printf("  <Y1[%i],Y1[%i]> = %g\n", i, j, vec2d_bDot(Y4[i],Y4[j]));
  printf("\n");
  printf("dimensions: (10, 1000000)\n");
  printf("testing time: %g\n", rtime);  
  
  vec2d_b** Y5 = (vec2d_b**) malloc(5 * sizeof(vec2d_b*));
  for (i=0; i<5; i++)
    Y5[i] = vec2d_bRandom(100000, 100);
  stime = get_time();
  if (GramSchmidt2d_b(Y5, 5))
    printf("GramSchmidt2d_b returned error\n");
  ftime = get_time();
  rtime = ftime-stime;
  printf("Resulting vectors should be orthonormal:\n");
  for (i=0; i<5; i++)
    for (j=i; j<5; j++)
      printf("  <Y1[%i],Y1[%i]> = %g\n", i, j, vec2d_bDot(Y5[i],Y5[j]));
  printf("\n");
  printf("dimensions: (100000, 100)\n");
  printf("testing time: %g\n", rtime);  
 
  vec2d_b** Y6 = (vec2d_b**) malloc(5 * sizeof(vec2d_b*));
  for (i=0; i<5; i++)
    Y6[i] = vec2d_bRandom(1000000, 10);
  stime = get_time();
  if (GramSchmidt2d_b(Y6, 5))
    printf("GramSchmidt2d_b returned error\n");
  ftime = get_time();
  rtime = ftime-stime;
  printf("Resulting vectors should be orthonormal:\n");
  for (i=0; i<5; i++)
    for (j=i; j<5; j++)
      printf("  <Y1[%i],Y1[%i]> = %g\n", i, j, vec2d_bDot(Y6[i],Y6[j]));
  printf("\n");
  printf("dimensions: (1000000, 10)\n");
  printf("testing time: %g\n", rtime);  

  /* clean up */
  vec2d_bDestroy(a);
  vec2d_bDestroy(b);
  vec2d_bDestroy(c);
  vec2d_bDestroy(d);
  vec2d_bDestroy(f);
  for (i=0; i<5; i++)
    vec2d_bDestroy(X[i]);
  for (i=0; i<5; i++)
    vec2d_bDestroy(Y1[i]);
  free(Y1);
  for (i=0; i<5; i++)
    vec2d_bDestroy(Y2[i]);
  free(Y2);
  for (i=0; i<5; i++)
    vec2d_bDestroy(Y3[i]);
  free(Y3);
  for (i=0; i<5; i++)
    vec2d_bDestroy(Y4[i]);
  free(Y4);
  for (i=0; i<5; i++)
    vec2d_bDestroy(Y5[i]);
  free(Y5);
  for (i=0; i<5; i++)
    vec2d_bDestroy(Y6[i]);
  free(Y6);
  
  return 0;
} /* end main */

