#include<math.h>
#include "f2c.h"

#define SIZE 10

extern int cbesj(complex*, real*, integer*, integer*, complex*, integer*,
		 integer*);

/* Using Recurrence formula of the Bessel function.  */
/* Jn'(x) = (Jn-1(x) - Jn+1(x))/2		     */

int cbesjp(complex *x, real *alpha, integer *kode, integer *n, complex *cy,
	  integer *nz, integer *ierr)
{
  complex cy1[SIZE], cy2[SIZE];
  real alpha1, alpha2;
  int i;

  alpha1 = *alpha - 1;
  alpha2 = *alpha + 1;

  cbesj(x, &alpha1, kode, n, cy1, nz, ierr);
  cbesj(x, &alpha2, kode, n, cy2, nz, ierr);

  for (i=0; i<*n; i++)
    {
      cy[i].r = (cy1[i].r - cy2[i].r)/2;
      cy[i].i = (cy1[i].i - cy2[i].i)/2;
    }

  return 0;
}
