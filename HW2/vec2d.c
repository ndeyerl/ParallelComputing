/* Nicole Deyerl
   SMU Mathematics
   Math 4370/6370
   3 March 2017 */ 

/* Inclusions */
#include "vec2d.h"


/* This file implements the operations defined in vec2d.h */

/* constructor (initializes values to 0.0) */
vec2d* vec2dNew(long int m, long int n) {

  long int i;
  
  /* if len is illegal return NULL pointer */
  if (m < 1 || n < 1) {
    fprintf(stderr, "vec2dNew error, illegal vector length = %li %li\n", m, n);
    return NULL;
  }
  
  /* otherwise try to allocate new vec2d object (return NULL on any failure) */
  vec2d *x = (vec2d *) malloc(sizeof(vec2d));
  if (x == NULL) return NULL;
  x->length1 = m;
  x->length2 = n;
  x->data = (double **) malloc( m*sizeof(double*) );
  for (i=0; i<m; i++)
    x->data[i] = (double *) calloc(n, sizeof(double));
  if (x->data == NULL) {
    free(x);
    return NULL;
  }
  return x;
};


/* destructor */
void vec2dDestroy(vec2d* v) {
  free(v->data);
  v->length1 = 0;
  v->length2 = 0;
}


/* write a vec2d to stdout */
int vec2dWrite(vec2d* v) {
  long int i;
  long int j;
  /* return with failure if data array isn't allocated */
  if (v->data == NULL) {
    fprintf(stderr, "vec2dWrite error, empty data array\n");
    return 1;
  }

  /* print data to screen and return */
  for (i=0; i<v->length1; i++){
	  for (j=0; j<v->length2; j++){
		  printf("  %.16g\n",v->data[i][j]);
	  }
  }
  printf("\n");
  return 0;
}


/* write a vec2d to a file */
int vec2dWriteFile(vec2d* v, const char* outfile) {
  long int i;
  long int j;
  FILE *fptr = NULL;

  /* return with failure if data array isn't allocated */
  if (v->data == NULL) {
    fprintf(stderr, "vec2dWriteFile error, empty data array\n");
    return 1;
  }

  /* return with failure if 'outfile' is empty */
  if (strlen(outfile) < 1) {
    fprintf(stderr, "vec2dWriteFile error, empty outfile\n");
    return 1;
  }

  /* open output file */
  fptr = fopen(outfile, "w");
  if (fptr == NULL) {
    fprintf(stderr, "vec2dWriteFile error, unable to open %s for writing\n", outfile);
    return 1;
  }

  /* print data to file */
  for (i=0; i<v->length1; i++){
	  for (j=0; j<v->length2; j++){
		  fprintf(fptr, "  %.16g\n",v->data[i][j]);
	  }
	  fprintf(fptr, "\n");
  }

  /* close output file and return */
  fclose(fptr);
  return 0;
}


/******** Arithmetic operations defined on a given vec2d ********/

/* x = a*y + b*z */
int vec2dLinearSum(vec2d* x, double a, vec2d* y, double b, vec2d* z) {
  long int i;
  long int j;

  /* check that array sizes match */
  if ((y->length1 != x->length1) || (z->length1 != x->length1) 
			|| (y->length2 != x->length2) || (z->length2 != x->length2)) {
    fprintf(stderr,"vec2dLinearSum error, vector sizes do not match\n");
    return 1;
  }
  
  /* check that data is not NULL */
  if (x->data == NULL || y->data == NULL || z->data == NULL) {
    fprintf(stderr, "vec2dLinearSum error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  x->data[i][j] = a*y->data[i][j] + b*z->data[i][j];
	  }
  }

  return 0;
}


/*   x = x*a  (scales x by a) */
int vec2dScale(vec2d* x, double a) {
  long int i;
  long int j;

  /* check that data is not NULL */
  if (x->data == NULL) {
    fprintf(stderr, "vec2dScale error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  x->data[i][j] *= a;
	  }
  }

  return 0;
}


/*   x = y  (copies y into x) */
int vec2dCopy(vec2d* x, vec2d* y) {
  long int i;
  long int j;
  
  /* check that array sizes match */
  if ((y->length1 != x->length1) || (y->length2 != x->length2)) {
    fprintf(stderr,"vec2dCopy error, vector sizes do not match\n");
    return 1;
  }
  
  /* check that data is not NULL */
  if (x->data == NULL || y->data == NULL) {
    fprintf(stderr, "vec2dCopy error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  x->data[i][j] = y->data[i][j];
	  }
  }

  return 0;
}


/*   x = a  (sets all entries of x to the scalar a) */
int vec2dConstant(vec2d* x, double a) {
  long int i;
  long int j;

  /* check that data is not NULL */
  if (x->data == NULL) {
    fprintf(stderr, "vec2dConst error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  x->data[i][j] = a;
	  }
  }
  return 0;
}






/******** scalar quantities derived from vectors ********/


/* min x_i */
double vec2dMin(vec2d* x) {
  double mn;
  long int i;
  long int j;
  mn = x->data[0][0];
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  mn = (mn < x->data[i][j]) ? mn : x->data[i][j];
	  }
  }
  return mn;
}


/* max x_i */
double vec2dMax(vec2d* x) {
  double mx;
  long int i;
  long int j;
  mx = x->data[0][0];
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  mx = (mx > x->data[i][j]) ? mx : x->data[i][j];
	  }
  }
  return mx;
}


/* dot-product of x and y */
double vec2dDot(vec2d* x, vec2d* y) {
  double sum;
  long int i;
  long int j;

  /* check that array sizes match */
  if ((y->length1 != x->length1) || (y->length2 != x->length2)){
    fprintf(stderr,"vec2dDot error, vector sizes do not match\n");
    return 0.0;
  }
  
  /* perform operation and return */
  sum = 0.0;
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  sum += (x->data[i][j] * y->data[i][j]);
	  }
  }
  return sum;
}


/* ||x||_2 */
double vec2dTwoNorm(vec2d* x) {
  double sum;
  long int i;
  long int j;
  sum = 0.0;
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  sum += (x->data[i][j] * x->data[i][j]);
	  }
  }
  return sqrt(sum);
}


/* ||x||_RMS */
double vec2dRmsNorm(vec2d* x) {
  double sum;
  long int i;
  long int j;
  sum = 0.0;
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  sum += (x->data[i][j] * x->data[i][j]);
	  }
  }
  return sqrt(sum/(x->length1 * x->length2));
}


/* ||x||_infty */
double vec2dMaxNorm(vec2d* x) {
  double mx;
  long int i;
  long int j;
  mx = 0.0;
  for (i=0; i<x->length1; i++){
	  for (j=0; j<x->length2; j++){
		  mx = (mx > fabs(x->data[i][j])) ? mx : fabs(x->data[i][j]);
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
vec2d* vec2dLinspace(double a, double b, long int m, long int n) {
  vec2d* x;
  long int i;
  long int j;

  /* create new vector of desired length */
  x = vec2dNew(m, n);
  if (x == NULL)  return NULL;

  /* fill in entries and return */
  /* works by filling each data set with the same linspace */
  for (i=0; i<m; i++){
	  for (j=0; j<n; j++){
		  x->data[i][j] = a + ((b-a)/(n*m-1))*j + i*n;
	  }
  }
  return x;
}


/* create a vector of uniformly-distributed random data */
vec2d* vec2dRandom(long int m, long int n) {
  vec2d* x;
  long int i;
  long int j;

  /* create new vector of desired length */
  x = vec2dNew(m, n);
  if (x == NULL)  return NULL;

  /* fill in entries and return */
  for (i=0; i<m; i++){
	  for (j=0; j<n; j++){
		  x->data[i][j] = random() / (pow(2.0,31.0) - 1.0);
	  }
  }
  return x;
}
