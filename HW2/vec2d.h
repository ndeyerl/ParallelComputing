/* Nicole Deyerl
   SMU Mathematics
   Math 4370/6370
   3 March 2017 */

#ifndef VEC2D_DEFINED__
#define VEC2D_DEFINED__

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


/* general utility macros */
#define True  1
#define False 0


/* This defines a simple arithmetic vector structure */
typedef struct _vec2d {
  long int length1; //m dim of m x n
  long int length2; //n dim of m x n
  double** data;
} vec2d;


/* General operations on vector structures */

/*   constructor (initializes values to 0.0) */
vec2d* vec2dNew(long int m, long int n);

/*   destructor */
void vec2dDestroy(vec2d* v);

/*   write a vec2d to stdout */
int vec2dWrite(vec2d* v);

/*   write a vec2d to a file */
int vec2dWriteFile(vec2d* v, const char* outfile);

/* arithmetic operations defined on a vec2d */
int vec2dLinearSum(vec2d* x, double a, vec2d* y,   /* x = a*y + b*z */
                   double b, vec2d* z);   
int vec2dScale(vec2d* x, double a);                /* x = x*a */
int vec2dCopy(vec2d* x, vec2d* y);                 /* x = y */
int vec2dConstant(vec2d* x, double a);             /* x = a */

/* scalar quantities derived from vectors */
double vec2dMin(vec2d* x);
double vec2dMax(vec2d* x);
double vec2dDot(vec2d* x, vec2d* y);
double vec2dTwoNorm(vec2d* x);
double vec2dRmsNorm(vec2d* x);
double vec2dMaxNorm(vec2d* x);

/* extra constructors */
vec2d* vec2dLinspace(double a, double b, long int m, long int n);
vec2d* vec2dRandom(long int m, long int n);

#endif
