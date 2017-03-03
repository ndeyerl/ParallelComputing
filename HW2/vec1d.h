/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   10 February 2017 */

#ifndef VEC1D_DEFINED__
#define VEC1D_DEFINED__

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


/* general utility macros */
#define True  1
#define False 0


/* This defines a simple arithmetic vector structure */
typedef struct _vec1d {
  long int length;
  double* data;
} vec1d;


/* General operations on vector structures */

/*   constructor (initializes values to 0.0) */
vec1d* vec1dNew(long int len);

/*   destructor */
void vec1dDestroy(vec1d* v);

/*   write a vec1d to stdout */
int vec1dWrite(vec1d* v);

/*   write a vec1d to a file */
int vec1dWriteFile(vec1d* v, const char* outfile);

/* arithmetic operations defined on a vec1d */
int vec1dLinearSum(vec1d* x, double a, vec1d* y,   /* x = a*y + b*z */
                   double b, vec1d* z);   
int vec1dScale(vec1d* x, double a);                /* x = x*a */
int vec1dCopy(vec1d* x, vec1d* y);                 /* x = y */
int vec1dConstant(vec1d* x, double a);             /* x = a */

/* scalar quantities derived from vectors */
double vec1dMin(vec1d* x);
double vec1dMax(vec1d* x);
double vec1dDot(vec1d* x, vec1d* y);
double vec1dTwoNorm(vec1d* x);
double vec1dRmsNorm(vec1d* x);
double vec1dMaxNorm(vec1d* x);

/* extra constructors */
vec1d* vec1dLinspace(double a, double b, long int n);
vec1d* vec1dRandom(long int n);

#endif
