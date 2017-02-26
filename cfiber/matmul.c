#include<stdio.h>
#include<math.h>
#include "f2c.h"

#define NCB  4
#define NRA  4
#define NCA  4

int matmul(real matA[][NCA], real matB[][NCB], real matC[][NRA])
{
  int i,j,k;

  for (k=0; k<NCB; k++)
    {
      for (i=0; i<NRA; i++)
	{
	  matC[i][k] = 0.0;
	  for (j=0; j<NCA; j++)
	    {
	      matC[i][k] = matC[i][k] + (matA[i][j] * matB[j][k]);
	    }
	}
    }
  return (0);
}
