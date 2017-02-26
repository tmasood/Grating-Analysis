#include<stdio.h>
#include<math.h>
#include "f2c.h"
 
#define PI M_PI
#define SIZE 4
#define DIMENSION 200
#define MAXSIZE 100

extern int besj(real*, real*, int*, real*, int*);
extern int besy(real*, real*, int*, real*);
extern int besk(real*, real*, int*, int*, real*, int*);
extern int besi(real*, real*, int*, int*, real*, int*);

/* find Coefficient transfer matrix for a given layer */
int fieldez(real betaroot, real k0, real nu, real *xF, int *wregion, real *Qtrans,
	    int nbound, real *typewave, real coefMat[][MAXSIZE], real *Ez)
{
  int ztype[DIMENSION];
  real alpha;
  real z2;
  real fnu;
  real besjqt[10], besyqt[10];
  real beskqt[10], besiqt[10];
  real Qt;
  
  int numx;
  int numr;
  int region;
  int i;
  int nz, n;
  int kode;

  n = 10;
  kode = 1;
  alpha = 1.0;
  fnu = 1.0;
  numx = DIMENSION;
  numr = nbound + 1;
  i = 0;
  while(i < numx)
    {
      region = wregion[i];
      Qt = Qtrans[region];
      z2 = Qt*xF[i];

      // kappa, J and Y
      if (typewave[region] > 0)
	{
	  besj(&z2, &alpha, &n, besjqt, &nz);
	  if (region == 0)
	    {
	      Ez[i] = coefMat[0][region]*besjqt[0];
	    }
	  else
	    {
	      besy(&z2, &fnu, &n, besyqt);
	      Ez[i] = coefMat[0][region]*besjqt[0] + coefMat[1][region]*besyqt[0];
	    }
	  ztype[i] = 19;
	}
      else
	{
	  while(1)
	    {
	      besi(&z2, &alpha, &kode, &n, besiqt, &nz);
	      if (region == 0)
		{
		  Ez[i] = coefMat[1][region]*besiqt[0];
		  break;
		}
	      besk(&z2, &fnu, &kode, &n, beskqt, &nz);
	      if (region == (numr - 1))
		{
		  Ez[i] = coefMat[0][region]*beskqt[0];
		  break;
		}
	      
	      Ez[i] = coefMat[0][region]*beskqt[0] + coefMat[1][region]*besiqt[0];
	      break;
	    }
	  ztype[i] = 8;
	}
      i++;
    }
  return (0);
}
