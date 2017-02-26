#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"
 
#define PI M_PI
#define SIZE 4
#define DIMENSION 200
#define MAXSIZE 100

extern int cbesj(complex*, real*, integer*, integer*, complex*,  integer*,
                 integer*);
extern int cbesy(complex*, real*, integer*, integer*, complex*,  integer*,
                 complex*, integer*);
extern int cbesk(complex*, real*, integer*, integer*, complex*,  integer*,
                 integer*);
extern int cbesi(complex*, real*, integer*, integer*, complex*,  integer*,
                 integer*);
extern int cbesjp(complex*, real*, integer*, integer*, complex*, integer*,
                  integer*);
extern int cbesyp(complex*, real*, integer*, integer*, complex*, integer*,
                  complex*, integer*);
extern int cbesip(complex*, real*, integer*, integer*, complex*, integer*,
                  integer*);
extern int cbeskp(complex*, real*, integer*, integer*, complex*, integer*,
                  integer*);
 
/* find Coefficient transfer matrix for a given layer */
int cfieldhr(complex betaroot, complex k0, real nu, real *xF,
	     int *wregion, complex *Qtrans,
	     int nbound, real *typewave, complex coefMat[][MAXSIZE],
	     complex *index, complex *Hr)
{
  real alpha;
  complex z2;
  complex besjqf[10], besyqf[10];
  complex beskqf[10], besiqf[10];
  complex besjpqf[10], besypqf[10];
  complex beskpqf[10], besipqf[10];
  complex a0;
  complex Qt;
  complex p1, p2, p;
  complex czero;
  complex cwrk[10];
  
  int numx;
  int numr;
  int region;
  int i;
  integer nz, n;
  integer kode;
  integer ierr;
  integer two;

  czero.r = 0.0;
  czero.i = 0.0;

  two = 2;
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

      pow_ci(&p1, &index[region], &two);
      p2 = cscalar_prod(k0, nu);
      a0 = c_prod(p2, p1); 

      z2 = cscalar_prod(Qt, xF[i]);
      /* kappa, J and Y */
      if (typewave[region] > 0)
	{
          cbesj(&z2, &alpha, &kode, &n, besjqf, &nz, &ierr);
          cbesjp(&z2, &alpha, &kode, &n, besjpqf, &nz, &ierr);
	  if (region == 0)
	    {
	      p1 = c_prod(betaroot, Qt);
	      p2 = c_prod(p1, coefMat[2][region]);
              p = c_prod(p2, besjpqf[0]);
	      Hr[i] = p;

              p1 = c_prod(coefMat[0][region], besjqf[0]);
              p2 = c_prod(a0, p1);
              p = cscalar_div(p2, xF[i]);
              p = c_add(p, Hr[i]);
	      Hr[i] = p;

              pow_ci(&p1, &Qt, &two);
              c_div(&p, &Hr[i], &p1);
	      p = cscalar_prod(p, (-1.0));
	      Hr[i] = p;
	    }
	  else
	    {
	      cbesy(&z2, &nu, &kode, &n, besyqf, &nz, cwrk, &ierr);
              cbesyp(&z2, &nu, &kode, &n, besypqf, &nz, cwrk, &ierr);

	      p1 = c_prod(coefMat[2][region], besjpqf[0]);
	      p2 = c_prod(coefMat[3][region], besypqf[0]);
	      p = c_add(p1, p2);
	      p1 = c_prod(betaroot, Qt);
	      p = c_prod(p, p1);
	      Hr[i] = p;

	      p1 = c_prod(coefMat[0][region], besjqf[0]);
	      p2 = c_prod(coefMat[1][region], besyqf[0]);
	      p = c_add(p1, p2);
	      p1 = c_prod(a0, p);
	      p = cscalar_div(p1, xF[i]);
	      p = c_add(p, Hr[i]);
	      Hr[i] = p;

	      pow_ci(&p1, &Qt, &two);
	      c_div(&p, &Hr[i], &p1);
	      p = cscalar_prod(p, (-1.0));
	      Hr[i] = p;
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
		  p1 = c_prod(coefMat[3][region], besipqf[0]);
		  p2 = c_prod(betaroot, Qt);
		  p = c_prod(p1, p2);
		  Hr[i] = p;

		  p1 = c_prod(coefMat[1][region], besiqf[0]);
		  p2 = c_prod(a0, p1);
		  p = cscalar_div(p2, xF[i]);
		  p = c_add(p, Hr[i]);
		  Hr[i] = p;

		  pow_ci(&p1, &Qt, &two);
		  c_div(&p, &Hr[i], &p1);
		  Hr[i] = p;
		  break;
		}
              cbesk(&z2, &nu, &kode, &n, beskqf, &nz, &ierr);
              cbeskp(&z2, &nu, &kode, &n, beskpqf, &nz, &ierr);
	      if (region == (numr - 1))
		{
		  p1 = c_prod(coefMat[2][region], beskpqf[0]);
		  p2 = c_prod(betaroot, Qt);
		  p = c_prod(p1, p2);
		  Hr[i] = p;

		  p1 = c_prod(coefMat[0][region], beskqf[0]);
		  p2 = c_prod(a0, p1);
		  p = cscalar_div(p2, xF[i]);
		  p = c_add(p, Hr[i]);
		  Hr[i] = p;

		  pow_ci(&p1, &Qt, &two);
		  c_div(&p, &Hr[i], &p1);
		  Hr[i] = p;
		  break;
		}
	      p1 = c_prod(coefMat[2][region], beskpqf[0]);
	      p2 = c_prod(coefMat[3][region], besipqf[0]);
	      p = c_add(p1, p2);
	      p1 = c_prod(betaroot, Qt);
	      p = c_prod(p, p1);
	      Hr[i] = p;

	      p1 = c_prod(coefMat[0][region], beskqf[0]);
	      p2 = c_prod(coefMat[1][region], besiqf[0]);
	      p = c_add(p1, p2);
	      p1 = c_prod(a0, p);
	      p = cscalar_div(p1, xF[i]);
	      p = c_add(p, Hr[i]);
	      Hr[i] = p;

	      pow_ci(&p1, &Qt, &two);
	      c_div(&p, &Hr[i], &p1);
	      Hr[i] = p;
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
	  pow_ci(&p2, &index[0], &two);
	  p = c_prod(k0, p2);
	  p = c_prod(coefMat[0][0], p);
	  p = c_add(p1, p);
	  p2 = cscalar_prod(Qtrans[0], 2);
	  c_div(&p1, &p, &p2);
	  p = cscalar_prod(p1, (-1.0));
	  Hr[0] = p;
	}
      else
	{
	  p1 = c_prod(coefMat[3][0], betaroot);
	  pow_ci(&p2, &index[0], &two);
	  p = c_prod(k0, p2);
	  p = c_prod(coefMat[1][0], p);
	  p = c_add(p1, p);
	  p2 = cscalar_prod(Qtrans[0], 2);
	  c_div(&p1, &p, &p2);
	  Hr[0] = p1;
	}
    }
  else
    {
      Hr[0] = czero;
    }
  return (0);
}
