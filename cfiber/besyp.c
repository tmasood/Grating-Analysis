#include<stdio.h>
#include<math.h>
#include "f2c.h"

#define SIZE 10

extern int besy(real*, real*, int*, real*);

/* Using Recurrence formula of the Bessel function.  */
/* Yn'(x) = (Yn-1(x) - Yn+1(x))/2                    */
 
int besyp(real *x, real *fnu, int *n, real *y)
{
  real y1[SIZE], y2[SIZE];
  real fnu1, fnu2;
  int i;
 
  fnu1 = *fnu - 1;
  fnu2 = *fnu + 1;

  besy(x, &fnu1, n, y1);
  besy(x, &fnu2, n, y2);

  for (i=0; i<*n; i++)
    {
      y[i] = (y1[i] - y2[i])/2;
    }
 
  return (0);
}
