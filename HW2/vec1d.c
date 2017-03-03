/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   10 February 2017 */

/* Inclusions */
#include "vec1d.h"


/* This file implements the operations defined in vec.h */

/* constructor (initializes values to 0.0) */
vec1d* vec1dNew(long int len) {

  /* if len is illegal return NULL pointer */
  if (len < 1) {
    fprintf(stderr, "vec1dNew error, illegal vector length = %li\n", len);
    return NULL;
  }
  
  /* otherwise try to allocate new vec1d object (return NULL on any failure) */
  vec1d *x = (vec1d *) malloc(sizeof(vec1d));
  if (x == NULL) return NULL;
  x->length = len;
  x->data = (double *) calloc(len, sizeof(double));
  if (x->data == NULL) {
    free(x);
    return NULL;
  }
  return x;
};


/* destructor */
void vec1dDestroy(vec1d* v) {
  free(v->data);
  v->length = 0;
}


/* write a vec1d to stdout */
int vec1dWrite(vec1d* v) {
  long int i;
  /* return with failure if data array isn't allocated */
  if (v->data == NULL) {
    fprintf(stderr, "vec1dWrite error, empty data array\n");
    return 1;
  }

  /* print data to screen and return */
  for (i=0; i<v->length; i++)
    printf("  %.16g\n",v->data[i]);
  printf("\n");
  return 0;
}


/* write a vec1d to a file */
int vec1dWriteFile(vec1d* v, const char* outfile) {
  long int i;
  FILE *fptr = NULL;

  /* return with failure if data array isn't allocated */
  if (v->data == NULL) {
    fprintf(stderr, "vec1dWriteFile error, empty data array\n");
    return 1;
  }

  /* return with failure if 'outfile' is empty */
  if (strlen(outfile) < 1) {
    fprintf(stderr, "vec1dWriteFile error, empty outfile\n");
    return 1;
  }

  /* open output file */
  fptr = fopen(outfile, "w");
  if (fptr == NULL) {
    fprintf(stderr, "vec1dWriteFile error, unable to open %s for writing\n", outfile);
    return 1;
  }

  /* print data to file */
  for (i=0; i<v->length; i++)
    fprintf(fptr, "  %.16g\n",v->data[i]);
  fprintf(fptr, "\n");

  /* close output file and return */
  fclose(fptr);
  return 0;
}


/******** Arithmetic operations defined on a given vec1d ********/

/* x = a*x + b*y */
int vec1dLinearSum(vec1d* x, double a, vec1d* y, double b, vec1d* z) {
  long int i;

  /* check that array sizes match */
  if ((y->length != x->length) || (z->length != x->length)) {
    fprintf(stderr,"vec1dLinearSum error, vector sizes do not match\n");
    return 1;
  }
  
  /* check that data is not NULL */
  if (x->data == NULL || y->data == NULL || z->data == NULL) {
    fprintf(stderr, "vec1dLinearSum error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length; i++)  
    x->data[i] = a*y->data[i] + b*z->data[i];
  return 0;
}


/*   x = x*a  (scales x by a) */
int vec1dScale(vec1d* x, double a) {
  long int i;

  /* check that data is not NULL */
  if (x->data == NULL) {
    fprintf(stderr, "vec1dScale error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length; i++)  
    x->data[i] *= a;
  return 0;
}


/*   x = y  (copies y into x) */
int vec1dCopy(vec1d* x, vec1d* y) {
  long int i;

  /* check that array sizes match */
  if (y->length != x->length) {
    fprintf(stderr,"vec1dCopy error, vector sizes do not match\n");
    return 1;
  }
  
  /* check that data is not NULL */
  if (x->data == NULL || y->data == NULL) {
    fprintf(stderr, "vec1dCopy error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length; i++)  
    x->data[i] = y->data[i];
  return 0;
}


/*   x = a  (sets all entries of x to the scalar a) */
int vec1dConstant(vec1d* x, double a) {
  long int i;

  /* check that data is not NULL */
  if (x->data == NULL) {
    fprintf(stderr, "vec1dConst error: empty data array\n");
    return 1;
  }

  /* perform operation and return */
  for (i=0; i<x->length; i++)  
    x->data[i] = a;
  return 0;
}






/******** scalar quantities derived from vectors ********/


/* min x_i */
double vec1dMin(vec1d* x) {
  double mn;
  long int i;
  mn = x->data[0];
  for (i=0; i<x->length; i++)  
    mn = (mn < x->data[i]) ? mn : x->data[i];
  return mn;
}


/* max x_i */
double vec1dMax(vec1d* x) {
  double mx;
  long int i;
  mx = x->data[0];
  for (i=0; i<x->length; i++)  
    mx = (mx > x->data[i]) ? mx : x->data[i];
  return mx;
}


/* dot-product of x and y */
double vec1dDot(vec1d* x, vec1d* y) {
  double sum;
  long int i;

  /* check that array sizes match */
  if (y->length != x->length) {
    fprintf(stderr,"vec1dDot error, vector sizes do not match\n");
    return 0.0;
  }
  
  /* perform operation and return */
  sum = 0.0;
  for (i=0; i<x->length; i++)  
    sum += (x->data[i] * y->data[i]);
  return sum;
}


/* ||x||_2 */
double vec1dTwoNorm(vec1d* x) {
  double sum;
  long int i;
  sum = 0.0;
  for (i=0; i<x->length; i++)  
    sum += (x->data[i] * x->data[i]);
  return sqrt(sum);
}


/* ||x||_RMS */
double vec1dRmsNorm(vec1d* x) {
  double sum;
  long int i;
  sum = 0.0;
  for (i=0; i<x->length; i++)  
    sum += (x->data[i] * x->data[i]);
  return sqrt(sum/x->length);
}


/* ||x||_infty */
double vec1dMaxNorm(vec1d* x) {
  double mx;
  long int i;
  mx = 0.0;
  for (i=0; i<x->length; i++)  
    mx = (mx > fabs(x->data[i])) ? mx : fabs(x->data[i]);
  return mx;
}


/******** extra constructors ********/

/* create a vector of linearly spaced data */
vec1d* vec1dLinspace(double a, double b, long int n) {
  vec1d* x;
  long int i;

  /* create new vector of desired length */
  x = vec1dNew(n);
  if (x == NULL)  return NULL;

  /* fill in entries and return */
  for (i=0; i<n; i++)
    x->data[i] = a + (b-a)/(n-1)*i;
  return x;
}


/* create a vector of uniformly-distributed random data */
vec1d* vec1dRandom(long int n) {
  vec1d* x;
  long int i;

  /* create new vector of desired length */
  x = vec1dNew(n);
  if (x == NULL)  return NULL;

  /* fill in entries and return */
  for (i=0; i<n; i++)
    x->data[i] = random() / (pow(2.0,31.0) - 1.0);
  return x;
}
