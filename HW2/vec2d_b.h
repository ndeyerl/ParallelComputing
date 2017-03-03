/* Nicole Deyerl
   SMU Mathematics
   Math 4370/6370
   3 March 2017 */

#ifndef VEC2D_B_DEFINED__
#define VEC2D_B_DEFINED__

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


/* general utility macros */
#define True  1
#define False 0


/* This defines a simple arithmetic vector structure */
typedef struct _vec2d_b {
  long int length1; //m dim of m x n
  long int length2; //n dim of m x n
  double* data;
} vec2d_b;


/* General operations on vector structures */

/*   constructor (initializes values to 0.0) */
vec2d_b* vec2d_bNew(long int m, long int n);

/*   destructor */
void vec2d_bDestroy(vec2d_b* v);

/*   write a vec2d to stdout */
int vec2d_bWrite(vec2d_b* v);

/*   write a vec2d to a file */
int vec2d_bWriteFile(vec2d_b* v, const char* outfile);

/* arithmetic operations defined on a vec2d */
int vec2d_bLinearSum(vec2d_b* x, double a, vec2d_b* y,   /* x = a*y + b*z */
                   double b, vec2d_b* z);   
int vec2d_bScale(vec2d_b* x, double a);                /* x = x*a */
int vec2d_bCopy(vec2d_b* x, vec2d_b* y);                 /* x = y */
int vec2d_bConstant(vec2d_b* x, double a);             /* x = a */

/* scalar quantities derived from vectors */
double vec2d_bMin(vec2d_b* x);
double vec2d_bMax(vec2d_b* x);
double vec2d_bDot(vec2d_b* x, vec2d_b* y);
double vec2d_bTwoNorm(vec2d_b* x);
double vec2d_bRmsNorm(vec2d_b* x);
double vec2d_bMaxNorm(vec2d_b* x);

/* extra constructors */
vec2d_b* vec2d_bLinspace(double a, double b, long int m, long int n);
vec2d_b* vec2d_bRandom(long int m, long int n);

#endif
