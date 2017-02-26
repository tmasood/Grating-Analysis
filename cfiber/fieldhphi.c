#include<stdio.h>
#include<math.h>
#include "f2c.h"

extern int besj(real*, real*, int*, real*, int*);
extern int besyp(real*, real*, int*, real*);
extern int besjp(real*, real*, int*, real*, int*);
extern int besy(real*, real*, int*, real*);
extern int besip(real*, real*, int*, int*, real*, int*);
extern int besk(real*, real*, int*, int*, real*, int*);
extern int besi(real*, real*, int*, int*, real*, int*);
extern int beskp(real*, real*, int*, int*, real*, int*);
 
#define PI M_PI
#define SIZE 4
#define DIMENSION 200
#define MAXSIZE 100
 
/* find Coefficient transfer matrix for a given layer */
int fieldhphi(real betaroot, real k0, real nu, real *xF, int *wregion, real *Qtrans,
              int nbound, real *typewave, real coefMat[][MAXSIZE],
	      real *index, real *Hphi)
{
  real alpha;
  real z2;
  real fnu;
  real besjqf[10], besyqf[10];
  real beskqf[10], besiqf[10];
  real besjpqf[10], besypqf[10];
  real beskpqf[10], besipqf[10];
  real a0, b0;
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
      a0 = (betaroot*nu)/((pow(Qt,2))*xF[i]);
      b0 = (pow(index[region],2)*k0)/Qt;
      z2 = Qt*xF[i];
      if (typewave[region] > 0)
	{
          besj(&z2, &alpha, &n, besjqf, &nz);
          besjp(&z2, &alpha, &n, besjpqf, &nz);  
	  if (region == 0)
	    {
	      Hphi[i] = coefMat[0][region]*b0*besjpqf[0] +
		coefMat[2][region]*a0*besjqf[0];
	    }
	  else
	    {
              besy(&z2, &fnu, &n, besyqf);
              besyp(&z2, &fnu, &n, besypqf);
	      Hphi[i] = coefMat[0][region]*b0*besjpqf[0] +
		coefMat[2][region]*a0*besjqf[0];
	      Hphi[i] += coefMat[1][region]*b0*besypqf[0] +
		coefMat[3][region]*a0*besyqf[0];
	    }
	}
      else
	{
	  /* kappa, K and I */
	  while(1)
	    {
              besi(&z2, &alpha, &kode, &n, besiqf, &nz);
              besip(&z2, &alpha, &kode, &n, besipqf, &nz);
	      if (region == 0)
		{
		  Hphi[i] = -(coefMat[1][region]*b0*besipqf[0]) -
		    coefMat[3][region]*a0*besiqf[0];
		  break;
		}
              besk(&z2, &fnu, &kode, &n, beskqf, &nz);
              beskp(&z2, &fnu, &kode, &n, beskpqf, &nz);
	      if (region == (numr - 1))
		{
		  Hphi[i] = -coefMat[0][region]*b0*beskpqf[0] -
		    coefMat[2][region]*a0*beskqf[0];
		  break;
		}
	      Hphi[i] = -coefMat[0][region]*b0*beskpqf[0] -
		coefMat[2][region]*a0*beskqf[0];
	      Hphi[i] -= coefMat[1][region]*b0*besipqf[0] +
		coefMat[3][region]*a0*besiqf[0];
	      break;
	    }
	}
      i++;
    }
  if (nu == 1)
    {
      if (typewave[0] > 0)
	{
	  Hphi[0] = coefMat[2][0]*betaroot + 
	    (coefMat[0][0]*k0*(pow(index[0],2)))/(2*Qtrans[0]);
	}
      else
	{
	  Hphi[0]=-(coefMat[3][0]*betaroot + 
	    (coefMat[1][0]*k0*(pow(index[0],2)))/(2*Qtrans[0]));
	}
    }
  else
    {
      Hphi[0] = 0;
    }
  return (0);
}
