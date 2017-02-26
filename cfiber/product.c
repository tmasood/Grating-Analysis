#include<stdio.h>
#include<math.h>
#include "f2c.h"

extern int transfer(real, real, real, real, real, real, real*);
extern int matmul(real*, real*, real*);

#define PI M_PI
#define SIZE 4

/* find total matrix product as function of beta */
int product(int nbound, real beta, real nu, real k0, real *index,
	    real *wrb, real M_product[][SIZE])
{
  int i, j, k;
  real prod[SIZE][SIZE];
  real ftm[SIZE][SIZE];

 /* initialize the prod and M_product matrices */
  for (i=0; i<SIZE; i++)
    {
      for (j=0; j<SIZE; j++)
	{
	  if (i == j)
	    {
	      M_product[i][j] = 1.0;
	      prod[i][j] = 1.0;
	    }
	  else
	    {
	      M_product[i][j] = 0.0;
	      prod[i][j] = 0.0;
	    }
	}
    }

  if (nbound == 1)
    {
      for (i=0; i<4; i++)
	{
	  for (j=0; j<4; j++)
	    {
	      prod[i][j] = M_product[i][j];
	    }
	}
    }
  else
    {
      k = 1;
      while (k < nbound)
	{
	  /* constructs 4x4 field transfer matrix for a given layer */
	  transfer(beta, nu, index[k], wrb[k], wrb[k+1], k0, ftm);
	  matmul(ftm, prod, M_product);
	  for (i=0; i<SIZE; i++)
	    {
	      for (j=0; j<SIZE; j++)
		{
		  prod[i][j] = M_product[i][j];
		}
	    }
	  k++;
	}
    }
  return (0);
}
