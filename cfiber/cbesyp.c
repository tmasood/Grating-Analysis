#include<stdio.h>
#include<math.h>
#include "f2c.h"

#define SIZE 10

extern int cbesy(complex*, real*, integer*, integer*, complex*,
		 integer*, complex*, integer*);

/* Using Recurrence formula of the Bessel function.  */
/* Yn'(x) = (Yn-1(x) - Yn+1(x))/2                    */
 
int cbesyp(complex *x, real *fnu, integer *kode, integer *n, complex *cy,
	   integer *nz, complex *cwrk, integer *ierr)
{
  complex cy1[SIZE], cy2[SIZE];
  real fnu1, fnu2;
  int i;
 
  fnu1 = *fnu - 1;
  fnu2 = *fnu + 1;

  cbesy(x, &fnu1, kode, n, cy1, nz, cwrk, ierr);
  cbesy(x, &fnu2, kode, n, cy2, nz, cwrk, ierr);

  for (i=0; i<*n; i++)
    {
      cy[i].r = (cy1[i].r - cy2[i].r)/2;
      cy[i].i = (cy1[i].i - cy2[i].i)/2;
    }
 
  return (0);
}
