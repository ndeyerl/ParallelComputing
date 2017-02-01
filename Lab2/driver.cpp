/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   13 January 2013
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// prototypes
double one_norm(int n, double* u);
void vector_sum(int n, double* u, double* v, double* w);
void vector_difference(int n, double* u, double* v, double* w);
void vector_product(int n, double* u, double* v, double* w);


// Description: Simple example program that calls multiple functions.
int main() {

  // local variables
  int i;
  double u[10], v[10], w[10], res;

  // fill in vectors u and v and output to screen
  for (i=0; i<10; i++) {
    u[i] = (double) i;
    v[i] = (double) (10-i);
  }
  printf("   u = ");
  for (i=0; i<10; i++)  printf("%g ",u[i]);
  printf("\n");
  printf("   v = ");
  for (i=0; i<10; i++)  printf("%g ",v[i]);
  printf("\n");

  // call the one-norm function and output result to screen
  res = one_norm(10,u);
  printf("   one norm of u = %g\n",res);
  
  // call the vector sum subroutine and output result to screen
  vector_sum(10,u,v,w);
  printf("   sum = ");
  for (i=0; i<10; i++)  printf("%g ",w[i]);
  printf("\n");
  
  // call the vector difference subroutine and output result to screen
  vector_difference(10,u,v,w);
  printf("   difference = ");
  for (i=0; i<10; i++)  printf("%g ",w[i]);
  printf("\n");
  
  // call the vector product subroutine and output result to screen
  vector_product(10,u,v,w);
  printf("   product = ");
  for (i=0; i<10; i++)  printf("%g ",w[i]);
  printf("\n");
  
  return 0;
}
