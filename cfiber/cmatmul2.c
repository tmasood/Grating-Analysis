#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"

#define NCA  4

int cmatmul2(complex matA[][NCA], complex *matB, complex *matC)
{
  int i,j;
  complex p1, p;
  complex czero;

  czero.r = 0.0;
  czero.i = 0.0;

  for (i=0; i<NCA; i++)
    {
      matC[i] = czero;
      for (j=0; j<NCA; j++)
	{
	  p1 = c_prod(matA[i][j], matB[j]);
	  p = c_add(matC[i], p1);
	  matC[i] = p;
	}
    }
  return (0);
}
