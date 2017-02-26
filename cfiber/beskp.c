#include<math.h>
#include "f2c.h"

#define SIZE 10

extern int besk(real*, real*, int*, int*, real*, int*);

/* Using Recurrence formula of the Bessel function.  */
/* In'(x) = (In-1(x) - In+1(x))/2		     */

int beskp(real *x, real *fnu, int *kode, int *n, real *y, int *nz)
{
  real y1[SIZE], y2[SIZE];
  real fnu1, fnu2;
  int i;

  fnu1 = *fnu - 1;
  fnu2 = *fnu + 1;

  besk(x, &fnu1, kode, n, y1, nz);
  besk(x, &fnu2, kode, n, y2, nz);

  for (i=0; i<*n; i++)
    {
      y[i] = -(y1[i] + y2[i])/2;
    }

  return (0);
}
