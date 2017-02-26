#include<stdio.h>
#include<math.h>
#include "f2c.h"

extern int cmplxtransfer(complex, real, complex, real, real, complex,
			 complex*);
extern int cmatmul(complex*, complex*, complex*);

#define PI M_PI
#define SIZE 4

/* find total matrix product as function of beta */
int cproduct(int nbound, complex beta, real nu, complex k0, complex *index,
	    real *wrb, complex M_product[][SIZE])
{
  int i, j, k;
  complex prod[SIZE][SIZE];
  complex ftm[SIZE][SIZE];
  complex czero;
  complex cone;
  
  czero.r = 0.0;
  czero.i = 0.0;

  cone.r = 1.0;
  cone.i = 0.0;

 /* initialize the prod and M_product matrices */
  for (i=0; i<SIZE; i++)
    {
      for (j=0; j<SIZE; j++)
	{
	  if (i == j)
	    {
	      M_product[i][j] = cone;
	      prod[i][j] = cone;
	    }
	  else
	    {
	      M_product[i][j] = czero;
	      prod[i][j] = czero;
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
	  cmplxtransfer(beta, nu, index[k], wrb[k], wrb[k+1], k0, ftm);
	  cmatmul(ftm, prod, M_product);
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
