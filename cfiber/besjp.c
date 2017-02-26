#include<math.h>
#include "f2c.h"

#define SIZE 10

extern int besj(real*, real*, int*, real*, int*);

/* Using Recurrence formula of the Bessel function.  */
/* Jn'(x) = (Jn-1(x) - Jn+1(x))/2		     */

int besjp(real *x, real *alpha, int *n, real *y, int *nz)
{
  real y1[SIZE], y2[SIZE];
  real alpha1, alpha2;
  int i;

  alpha1 = *alpha - 1;
  alpha2 = *alpha + 1;

  besj(x, &alpha1, n, y1, nz);
  besj(x, &alpha2, n, y2, nz);

  for (i=0; i<*n; i++)
    {
      y[i] = (y1[i] - y2[i])/2;
    }

  return 0;
}
