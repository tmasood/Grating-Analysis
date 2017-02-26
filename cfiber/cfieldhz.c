#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"
 
#define PI M_PI
#define SIZE 4
#define DIMENSION 200
#define MAXSIZE 100

extern int cbesj(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbesy(complex*, real*, integer*, integer*, complex*, integer*,
		 complex*, integer*);
extern int cbesk(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbesi(complex*, real*, integer*, integer*, complex*, integer*, integer*);

/* find Coefficient transfer matrix for a given layer */
int cfieldhz(complex betaroot, complex k0, real nu, real *xF, int *wregion, complex *Qtrans,
	    int nbound, real *typewave, complex coefMat[][MAXSIZE], complex *Hz)
{
  real alpha;
  complex z2;
  real fnu;
  complex besjqt[10], besyqt[10];
  complex beskqt[10], besiqt[10];
  complex Qt;

  complex cwrk[10];
  complex p1, p2, p;

  int numx;
  int numr;
  int region;
  int i;
  integer nz, n;
  integer kode;
  integer ierr;

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
      z2 = cscalar_prod(Qt, xF[i]);
      
      // kappa, J and Y
      if (typewave[region] > 0)
	{
	  cbesj(&z2, &alpha, &kode, &n, besjqt, &nz, &ierr);
	  if (region == 0)
	    {
	      Hz[i] = c_prod(coefMat[2][region], besjqt[0]);
	    }
	  else
	    {
	      cbesy(&z2, &fnu, &kode, &n, besyqt, &nz, cwrk, &ierr);
	      p1 = c_prod(coefMat[2][region], besjqt[0]);
	      p2 = c_prod(coefMat[3][region], besyqt[0]);
	      p = c_add(p1, p2);
	      Hz[i] = p;
	    }
	}
      else
	{
	  while(1)
	    {
	      cbesi(&z2, &alpha, &kode, &n, besiqt, &nz, &ierr);
	      if (region == 0)
		{
		  Hz[i] = c_prod(coefMat[3][region], besiqt[0]);
		  break;
		}
	      cbesk(&z2, &fnu, &kode, &n, beskqt, &nz, &ierr);
	      if (region == (numr - 1))
		{
		  Hz[i] = c_prod(coefMat[2][region], beskqt[0]);
		  break;
		}
	      p1 = c_prod(coefMat[2][region], beskqt[0]);
	      p2 = c_prod(coefMat[3][region], besiqt[0]);
	      p = c_add(p1, p2);
	      Hz[i] = p;
	      break;
	    }
	}
      i++;
    }
  return (0);
}
