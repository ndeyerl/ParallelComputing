/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   13 January 2013
*/

#include <math.h>

// Description: computes the 1-norm of the vector u (of length n)
double one_norm(int n, double* u) {

  double sum=0.0;
  for (int i=0; i<n; i++)
    sum += fabs(u[i]);

  return sum;
}
