/* Nicole Deyerl
   SMU Mathematics
   Math 4370/6370
   3 March 2017 */ 

/* Inclusions */
#include "vec2d_b.h"


/* This file implements the operations defined in vec2d_b.h */

/* constructor (initializes values to 0.0) */
vec2d_b* vec2d_bNew(long int m, long int n) {

  long int i;
  
  /* if len is illegal return NULL pointer */
  if (m < 1 || n < 1) {
    fprintf(stderr, "vec2d_bNew error, illegal vector length = %li %li\n", m, n);
    return NULL;
  }
  
  /* otherwise try to allocate new vec2d_b object (return NULL on any failure) */
  vec2d_b *x = (vec2d_b *) malloc(sizeof(vec2d_b));
  if (x == NULL) return NULL;
  x->length1 = m;
  x->length2 = n;
  x->data = (double *) calloc(m*n, sizeof(double));
  if (x->data == NULL) {
    free(x);
    return NULL;
  }
  return x;
};


/* destructor */
void vec2d_bDestroy(vec2d_b* v) {
  free(v->data);
  v->length1 = 0;
  v->length2 = 0;
}


/* write a vec2d_b to stdout */
int vec2d_bWrite(vec2d_b* v) {
  long int i;
  long int j;
  /* return with failure if data array isn't allocated */
  if (v->data == NULL) {
    fprintf(stderr, "vec2d_bWrite error, empty data array\n");
    return 1;
  }

  /* print data to screen and return */
  for (i=0; i<v->length1; i++){
	  for (j=0; j<v->length2; j++){
		  printf("  %.16g\n",v->data[i*v->length2 + j]);
	  }
  }
  printf("\n");
  return 0;
}


/* write a vec2d_b to a file */
int vec2d_bWriteFile(vec2d_b* v, const char* outfile) {
  long int i;
  long int j;
  FILE *fptr = NULL;

  /* return with failure if data array isn't allocated */
  if (v->data == NULL) {
    fprintf(stderr, "vec2d_bWriteFile error, empty data array\n");
    return 1;
  }

  /* return with failure if 'outfile' is empty */
  if (strlen(outfile) < 1) {
    fprintf(stderr, "vec2d_bWriteFile error, empty outfile\n");
    return 1;
  }

  /* open output file */
  fptr = fopen(outfile, "w");
  if (fptr == NULL) {
    fprintf(stderr, "vec2d_bWriteFile error, unable to open %s for writing\n", outfile);
    return 1;
  }

  /* print data to file */
  for (i=0; i<v->length1; i++){
	  for (j=0; j<v->length2; j++){
		  fprintf(fptr, "  %.16g\n",v->data[i*v->length2 + j]);
	  }
	  fprintf(fptr, "\n");
  }

  /* close output file and return */
  fclose(fptr);
  return 0;
}


/******** Arithmetic operations defined on a given vec2d_b ********/

/* x = a*y + b*z */
int vec2d_bLinearSum(vec2d_b* x, double a, vec2d_b* y, double b, vec2d_b* z) {
  long int i;
  long int j;

  /* check that array sizes match */
  if ((y->length1 != x->length1) || (z->length1 != x->length1) 
			|| (y->length2 != x->length2) || (z->length2 != x->length2)) {
    fprintf(stderr,"vec2d_bLinearSum error, vector sizes do not match\n");
    return 1;
  }
  
  /* check that data is not NULL */
  if (x->data == NULL || y->data == NULL || z->data == NULL) {
    fprintf(stderr, "vec2d_bLinearSum error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  x->data[i*x->length2 + j] = a*y->data[i*y->length2 + j] + b*z->data[i*x->length2 + j];
	  }
  }

  return 0;
}


/*   x = x*a  (scales x by a) */
int vec2d_bScale(vec2d_b* x, double a) {
  long int i;
  long int j;

  /* check that data is not NULL */
  if (x->data == NULL) {
    fprintf(stderr, "vec2d_bScale error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  x->data[i*x->length2 + j] *= a;
	  }
  }

  return 0;
}


/*   x = y  (copies y into x) */
int vec2d_bCopy(vec2d_b* x, vec2d_b* y) {
  long int i;
  long int j;
  
  /* check that array sizes match */
  if ((y->length1 != x->length1) || (y->length2 != x->length2)) {
    fprintf(stderr,"vec2d_bCopy error, vector sizes do not match\n");
    return 1;
  }
  
  /* check that data is not NULL */
  if (x->data == NULL || y->data == NULL) {
    fprintf(stderr, "vec2d_bCopy error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  x->data[i*x->length2 + j] = y->data[i*y->length2 + j];
	  }
  }

  return 0;
}


/*   x = a  (sets all entries of x to the scalar a) */
int vec2d_bConstant(vec2d_b* x, double a) {
  long int i;
  long int j;

  /* check that data is not NULL */
  if (x->data == NULL) {
    fprintf(stderr, "vec2d_bConst error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  x->data[i*x->length2 + j] = a;
	  }
  }
  return 0;
}






/******** scalar quantities derived from vectors ********/


/* min x_i */
double vec2d_bMin(vec2d_b* x) {
  double mn;
  long int i;
  long int j;
  mn = x->data[0];
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  mn = (mn < x->data[i*x->length2 + j]) ? mn : x->data[i*x->length2 + j];
	  }
  }
  return mn;
}


/* max x_i */
double vec2d_bMax(vec2d_b* x) {
  double mx;
  long int i;
  long int j;
  mx = x->data[0];
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  mx = (mx > x->data[i*x->length2 + j]) ? mx : x->data[i*x->length2 + j];
	  }
  }
  return mx;
}


/* dot-product of x and y */
double vec2d_bDot(vec2d_b* x, vec2d_b* y) {
  double sum;
  long int i;
  long int j;

  /* check that array sizes match */
  if ((y->length1 != x->length1) || (y->length2 != x->length2)){
    fprintf(stderr,"vec2d_bDot error, vector sizes do not match\n");
    return 0.0;
  }
  
  /* perform operation and return */
  sum = 0.0;
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  sum += (x->data[i*x->length2 + j] * y->data[i*y->length2 + j]);
	  }
  }
  return sum;
}


/* ||x||_2 */
double vec2d_bTwoNorm(vec2d_b* x) {
  double sum;
  long int i;
  long int j;
  sum = 0.0;
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  sum += (x->data[i*x->length2 + j] * x->data[i*x->length2 + j]);
	  }
  }
  return sqrt(sum);
}


/* ||x||_RMS */
double vec2d_bRmsNorm(vec2d_b* x) {
  double sum;
  long int i;
  long int j;
  sum = 0.0;
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  sum += (x->data[i*x->length2 + j] * x->data[i*x->length2 + j]);
	  }
  }
  return sqrt(sum/(x->length1 * x->length2));
}


/* ||x||_infty */
double vec2d_bMaxNorm(vec2d_b* x) {
  double mx;
  long int i;
  long int j;
  mx = 0.0;
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  mx = (mx > fabs(x->data[i*x->length2 + j])) ? mx : fabs(x->data[i*x->length2 + j]);
	  }
  }
  return mx;
}


/******** extra constructors ********/

/* create a vector of linearly spaced data */
/* fills vector row-wise (ie [0,8]-> 3x3 -> [0 1 2
 *                                           3 4 5
 *                                           6 7 8]
 */
vec2d_b* vec2d_bLinspace(double a, double b, long int m, long int n) {
  vec2d_b* x;
  long int i;
  long int j;

  /* create new vector of desired length */
  x = vec2d_bNew(m, n);
  if (x == NULL)  return NULL;

  /* fill in entries and return */
  /* works by filling each data set with the same linspace */
  for (i=0; i<m; i++){
	  for (j=0; j<n; j++){
		  x->data[i*x->length2 + j] = a + ((b-a)/(n*m-1))*j + i*n;
	  }
  }
  return x;
}


/* create a vector of uniformly-distributed random data */
vec2d_b* vec2d_bRandom(long int m, long int n) {
  vec2d_b* x;
  long int i;
  long int j;

  /* create new vector of desired length */
  x = vec2d_bNew(m, n);
  if (x == NULL)  return NULL;

  /* fill in entries and return */
  for (i=0; i<m; i++){
	  for (j=0; j<n; j++){
		  x->data[i*x->length2 + j] = random() / (pow(2.0,31.0) - 1.0);
	  }
  }
  return x;
}
