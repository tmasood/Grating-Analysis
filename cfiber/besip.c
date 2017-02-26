#include<math.h>
#include "f2c.h"

#define SIZE 10

extern int besi(real*, real*, int*, int*, real*, int*);

/* Using Recurrence formula of the Bessel function.  */
/* In'(x) = (In-1(x) - In+1(x))/2		     */

int besip(real *x, real *alpha, int *kode, int *n, real *y, int *nz)
{
  real y1[SIZE], y2[SIZE];
  real alpha1, alpha2;
  int i;

  alpha1 = *alpha - 1;
  alpha2 = *alpha + 1;

  besi(x, &alpha1, kode, n, y1, nz);
  besi(x, &alpha2, kode, n, y2, nz);

  for (i=0; i<*n; i++)
    {
      y[i] = (y1[i] + y2[i])/2;
    }

  return 0;
}
