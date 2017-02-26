#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"

#define NCB  4
#define NRA  4
#define NCA  4

int cmatmul(complex matA[][NCA], complex matB[][NCB], complex matC[][NRA])
{
  int i,j,k;

  complex p1, p2, p;
  complex czero;

  czero.r = 0.0;
  czero.i = 0.0;

  for (k=0; k<NCB; k++)
    {
      for (i=0; i<NRA; i++)
	{
	  matC[i][k] = czero;
	  for (j=0; j<NCA; j++)
	    {
	      p1 = c_prod(matA[i][j], matB[j][k]);
	      p = c_add(p1, matC[i][k]);
	      matC[i][k] = p;
	    }
	}
    }
  return (0);
}
