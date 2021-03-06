#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"

extern int cbesj(complex*, real*, integer*, integer*, complex*,  integer*,
                 integer*);
extern int cbesyp(complex*, real*, integer*, integer*, complex*, integer*,
                  complex*, integer*);
extern int cbesjp(complex*, real*, integer*, integer*, complex*, integer*,
                  integer*);
extern int cbesy(complex*, real*, integer*, integer*, complex*,  integer*,
                 complex*, integer*);
extern int cbesip(complex*, real*, integer*, integer*, complex*, integer*,
                  integer*);
extern int cbesk(complex*, real*, integer*, integer*, complex*,  integer*,
                 integer*);
extern int cbesi(complex*, real*, integer*, integer*, complex*,  integer*,
                 integer*);
extern int cbeskp(complex*, real*, integer*, integer*, complex*, integer*,
                  integer*);

#define PI M_PI
#define SIZE 4
#define DIMENSION 200
#define MAXSIZE 100
 
/* find Coefficient transfer matrix for a given layer */
int cfieldhphi(complex betaroot, complex k0, real nu, real *xF, int *wregion,
	      complex *Qtrans, int nbound, real *typewave,
	      complex coefMat[][MAXSIZE], complex *index, complex *Hphi)
{
  real alpha;
  complex z2;
  real fnu;
  complex besjqf[10], besyqf[10];
  complex beskqf[10], besiqf[10];
  complex besjpqf[10], besypqf[10];
  complex beskpqf[10], besipqf[10];
  complex a0, b0;
  complex Qt;
  complex p1, p2, p3, p;
  complex czero;
  complex cwrk[10];

  int numx;
  int numr;
  int region;
  int i;
  integer nz, n;
  integer kode, ierr;
  integer two;
 
  czero.r = 0.0;
  czero.i = 0.0;
  n = 10;
  two = 2;
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
      p1 = cscalar_prod(betaroot, nu);
      pow_ci(&p2, &Qt, &two);
      p = cscalar_prod(p2, xF[i]);
      c_div(&a0, &p1, &p);
      
      pow_ci(&p1, &index[region], &two);
      p2 = c_prod(p1, k0);
      c_div(&b0, &p2, &Qt);

      z2 = cscalar_prod(Qt, xF[i]);
      if (typewave[region] > 0)
	{
          cbesj(&z2, &alpha, &kode, &n, besjqf, &nz, &ierr);
          cbesjp(&z2, &alpha, &kode, &n, besjpqf, &nz, &ierr);  
	  if (region == 0)
	    {
	      p1 = c_prod(coefMat[0][region], b0);
	      p1 = c_prod(p1, besjpqf[0]);
	      p2 = c_prod(coefMat[2][region], a0);
	      p2 = c_prod(p2, besjqf[0]);
	      p = c_add(p1, p2);
	      Hphi[i] = p;
	    }
	  else
	    {
              cbesy(&z2, &fnu, &kode, &n, besyqf, &nz, cwrk, &ierr);
              cbesyp(&z2, &fnu, &kode, &n, besypqf, &nz, cwrk, &ierr);

	      p1 = c_prod(coefMat[0][region], b0);
	      p1 = c_prod(p1, besjpqf[0]);
	      p2 = c_prod(coefMat[2][region],a0);
	      p2 = c_prod(p2, besjqf[0]);
	      p = c_add(p1, p2);
	      Hphi[i] = p;
	      
	      p1 = c_prod(coefMat[1][region], b0);
	      p1 = c_prod(p1, besypqf[0]);
	      p2 = c_prod(coefMat[3][region], a0);
	      p2 = c_prod(p2, besyqf[0]);
	      p = c_add(p1, p2);
	      Hphi[i] = c_add(Hphi[i], p);
	    }
	}
      else
	{
	  /* kappa, K and I */
	  while(1)
	    {
              cbesi(&z2, &alpha, &kode, &n, besiqf, &nz, &ierr);
              cbesip(&z2, &alpha, &kode, &n, besipqf, &nz, &ierr);
	      if (region == 0)
		{
		  p1 = c_prod(coefMat[1][region], b0);
		  p1 = c_prod(p1, besipqf[0]);
		  p1 = cscalar_prod(p1, (-1.0));
		  p2 = c_prod(coefMat[3][region], a0);
		  p2 = c_prod(p2, besiqf[0]);
		  p = c_sub(p1, p2);
		  Hphi[i] = p;
		  break;
		}
              cbesk(&z2, &fnu, &kode, &n, beskqf, &nz, &ierr);
              cbeskp(&z2, &fnu, &kode, &n, beskpqf, &nz, &ierr);
	      if (region == (numr - 1))
		{
		  p1 = c_prod(coefMat[0][region], b0);
		  p1 = c_prod(p1, beskpqf[0]);
		  p1 = cscalar_prod(p1, (-1.0));
		  p2 = c_prod(coefMat[2][region], a0);
		  p2 = c_prod(p2, beskqf[0]);
		  p = c_sub(p1, p2);
		  Hphi[i] = p;
		  break;
		}
	      p1 = c_prod(coefMat[0][region], b0);
	      p1 = c_prod(p1, beskpqf[0]);
	      p1 = cscalar_prod(p1, (-1.0));
	      p2 = c_prod(coefMat[2][region], a0);
	      p2 = c_prod(p2, beskqf[0]);
	      p = c_prod(p1, p2);
	      Hphi[i] = p;

	      p1 = c_prod(coefMat[1][region], b0);
	      p1 = c_prod(p1, besipqf[0]);
	      p2 = c_prod(coefMat[3][region], a0);
	      p2 = c_prod(p2, besiqf[0]);
	      p = c_add(p1, p2);
	      Hphi[i] = c_sub(Hphi[i], p);
	      break;
	    }
	}
      i++;
    }
  if (nu == 1)
    {
      if (typewave[0] > 0)
	{
	  p1 = c_prod(coefMat[2][0], betaroot);
	  p2 = c_prod(coefMat[0][0], k0);
	  pow_ci(&p, &index[0], &two);
	  p2 = c_prod(p2, p);
	  p = cscalar_prod(Qtrans[0], 2);
	  c_div(&p3, &p2, &p);
	  Hphi[0] = c_add(p1, p3);
	}
      else
	{
	  p1 = c_prod(coefMat[3][0], betaroot);
	  p2 = c_prod(coefMat[1][0], k0);
	  pow_ci(&p, &index[0], &two);
	  p2 = c_prod(p2, p);
	  p = cscalar_prod(Qtrans[0], 2);
	  c_div(&p3, &p2, &p);
	  p = c_add(p3, p1);
	  p = cscalar_prod(p, (-1.0));
	  Hphi[0] = p;
	}
    }
  else
    {
      Hphi[0] = czero;
    }
  return (0);
}
