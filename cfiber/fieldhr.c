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
int fieldhr(real betaroot, real k0, real nu, real *xF, int *wregion, real *Qtrans,
	    int nbound, real *typewave, real coefMat[][MAXSIZE],
	    real *index, real *Hr)
{
  real alpha;
  real z2;
  real besjqf[10], besyqf[10];
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
      a0 = k0*nu*(pow(index[region],2));
      z2 = Qt*xF[i];      
      /* kappa, J and Y */
      if (typewave[region] > 0)
	{
          besj(&z2, &alpha, &n, besjqf, &nz);
          besjp(&z2, &alpha, &n, besjpqf, &nz);
	  if (region == 0)
	    {
	      Hr[i] = betaroot*Qt*(coefMat[2][region]*besjpqf[0]);
	      Hr[i]+= (a0*(coefMat[0][region]*besjqf[0]))/xF[i];
	      Hr[i] = -Hr[i]/pow(Qt,2);
	    }
	  else
	    {
	      besy(&z2, &nu, &n, besyqf);
              besyp(&z2, &nu, &n, besypqf);
	      Hr[i] = betaroot*Qt*((coefMat[2][region]*besjpqf[0]) +
		(coefMat[3][region]*besypqf[0]));
	      Hr[i] += (a0*((coefMat[0][region]*besjqf[0]) +
		(coefMat[1][region]*besyqf[0])))/xF[i];
	      Hr[i] = -Hr[i]/pow(Qt,2);
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
		  Hr[i] = betaroot*Qt*(coefMat[3][region]*besipqf[0]);
		  Hr[i]+= (a0*(coefMat[1][region]*besiqf[0]))/xF[i];
		  Hr[i] = Hr[i]/pow(Qt,2);
		  break;
		}
              besk(&z2, &nu, &kode, &n, beskqf, &nz);
              beskp(&z2, &nu, &kode, &n, beskpqf, &nz);
	      if (region == (numr - 1))
		{
		  Hr[i] = betaroot*Qt*(coefMat[2][region]*beskpqf[0]);
		  Hr[i]+= (a0*(coefMat[0][region]*beskqf[0]))/xF[i];
		  Hr[i] = Hr[i]/pow(Qt,2);
		  break;
		}
	      Hr[i] = betaroot*Qt*((coefMat[2][region]*beskpqf[0]) +
		(coefMat[3][region]*besipqf[0]));
	      Hr[i] += (a0*((coefMat[0][region]*beskqf[0]) +
			   (coefMat[1][region]*besiqf[0])))/xF[i];
	      Hr[i] = Hr[i]/pow(Qt,2);
	      break;
	    }
	}
      i++;
    }
  if (nu == 1)
    {
      if (typewave[0] > 0)
	{
	  Hr[0] = -((coefMat[2][0]*betaroot) + 
		    (coefMat[0][0]*k0*pow(index[0],2)))/(2*Qtrans[0]);
	}
      else
	{
	  Hr[0] = ((coefMat[3][0]*betaroot) + 
		   (coefMat[1][0]*k0*pow(index[0],2)))/(2*Qtrans[0]);
	}
    }
  else
    {
      Hr[0] = 0;
    }
  return (0);
}
