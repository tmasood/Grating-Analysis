#include<stdio.h>
#include<math.h>
#include "f2c.h"

#define NCA  4

int matmul2(real matA[][NCA], real *matB, real *matC)
{
  int i,j;

  for (i=0; i<NCA; i++)
    {
      matC[i] = 0.0;
      for (j=0; j<NCA; j++)
	{
	  matC[i] = matC[i] + (matA[i][j] * matB[j]);
	}
    }
  return (0);
}
