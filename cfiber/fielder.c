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
extern int besjp(real*, real*, int*, real*, int*);
extern int besyp(real*, real*, int*, real*);
extern int besip(real*, real*, int*, int*, real*, int*);
extern int beskp(real*, real*, int*, int*, real*, int*);
 
/* find Coefficient transfer matrix for a given layer */
int fielder(real betaroot, real k0, real nu, real *xF, int *wregion, real *Qtrans,
            int nbound, real *typewave, real coefMat[][MAXSIZE], real *Er)
{
  real alpha;
  real z2;
  real besjqf[10],besyqf[10];
  real beskqf[10], besiqf[10];
  real besjpqf[10], besypqf[10];
  real beskpqf[10], besipqf[10];
  real a0;
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
  numx = DIMENSION;
  numr = nbound + 1;
  i = 0;
  while(i < numx)
    {
      region = wregion[i];
      Qt = Qtrans[region];
      a0 = k0*nu;
      z2 = Qt*xF[i];      
      /* kappa, J and Y */
      if (typewave[region] > 0)
	{
          besj(&z2, &alpha, &n, besjqf, &nz);
          besjp(&z2, &alpha, &n, besjpqf, &nz);
	  if (region == 0)
	    {
	      Er[i] = betaroot*Qt*(coefMat[0][region]*besjpqf[0]);
	      Er[i] += (a0*(coefMat[2][region]*besjqf[0]))/xF[i];
	      Er[i] = Er[i]/pow(Qt,2);
	    }
	  else
	    {
              besyp(&z2, &nu, &n, besypqf);
              besy(&z2, &nu, &n, besyqf);
	      Er[i] = betaroot*Qt*((coefMat[0][region]*besjpqf[0]) +
		(coefMat[1][region]*besypqf[0]));
	      Er[i] += (a0*(coefMat[2][region]*besjqf[0] +
			   coefMat[3][region]*besyqf[0]))/xF[i];
	      Er[i] = Er[i]/pow(Qt,2);
	    }
	}
      else
	{
	  while(1)
	    {
              besi(&z2, &alpha, &kode, &n, besiqf, &nz);
              besip(&z2, &alpha, &kode, &n, besipqf, &nz);
	      if (region == 0)
		{
		  Er[i] = betaroot*Qt*(coefMat[1][region]*besipqf[0]);
		  Er[i]+= (a0*(coefMat[3][region]*besiqf[0]))/xF[i];
		  Er[i] = -Er[i]/pow(Qt,2);
		  break;
		}
              besk(&z2, &nu, &kode, &n, beskqf, &nz);
              beskp(&z2, &nu, &kode, &n, beskpqf, &nz);
	      if (region == (numr - 1))
		{
		  Er[i] = betaroot*Qt*(coefMat[0][region]*beskpqf[0]);
		  Er[i]+= (a0*(coefMat[2][region]*beskqf[0]))/xF[i];
		  Er[i] = -Er[i]/pow(Qt,2);
		  break;
		}
	      Er[i] = betaroot*Qt*((coefMat[0][region]*beskpqf[0]) +
		(coefMat[1][region]*besipqf[0]));
	      Er[i] += (a0*(coefMat[2][region]*beskqf[0] +
		coefMat[3][region]*besiqf[0]))/xF[i];
	      Er[i] = -Er[i]/pow(Qt,2);
	      break;
	    }
	}
      i++;
    }

  if (nu == 1.0)
    {
      if (typewave[0] > 0)
	{
	  Er[0] = (coefMat[0][0]*betaroot + coefMat[2][0]*k0)/(2*Qtrans[0]);
	}
      else
	{
	  Er[0] = -(coefMat[1][0]*betaroot + coefMat[3][0]*k0)/(2*Qtrans[0]);
	}
    }
  else
    {
      Er[0] = 0;
    }
  return (0);
}
