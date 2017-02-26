#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"

extern int cbesj(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbesyp(complex*, real*, integer*, integer*, complex*, integer*, complex*,
		  integer*);
extern int cbesjp(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbesy(complex*, real*, integer*, integer*, complex*, integer*, complex*,
		 integer*);
extern int cbesip(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbesk(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbesi(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbeskp(complex*, real*, integer*, integer*, complex*, integer*, integer*);

#define PI M_PI
#define SIZE 4

int cmplxtransfer(complex beta, real nu, complex nref, real r1, real r2,
	     complex k0, complex ftm[][SIZE])
{
  complex lay[SIZE][SIZE];
  complex kappa, gamma;
  complex alpha1;
  complex alpha2;
  complex rho;
  complex sigma;
  complex z1;
  complex z2;
  complex del1, del2, del3, del4;
  real alpha;
  real fnu;
  complex del1bjy[10], del1bjpy[10];
  complex del1byy[10], del1bypy[10];
  complex del2bjy[10], del2bjpy[10];
  complex del2byy[10], del2bypy[10];
  complex del3bjy1[10], del3bjy2[10];
  complex del3byy1[10], del3byy2[10];
  complex del4bjpy1[10], del4bjpy2[10];
  complex del4bypy1[10], del4bypy2[10];

  complex del1bipy[10], del1bky[10], del1biy[10];
  complex del1bkpy[10], del2bipy[10], del2bky[10],del2bkpy[10];
  complex del2biy[10], del3biy2[10], del3bky1[10];
  complex del3biy1[10], del3bky2[10], del4bipy2[10];
  complex del4bkpy1[10], del4bipy1[10], del4bkpy2[10];

  complex cwrk[10];

  complex p1, p2, p;
  complex czero;

  integer nz, n;
  integer kode;
  integer ierr;
  integer two;
  int i, j;
  int state;

  czero.r = 0.0;
  czero.i = 0.0;
  kode = 1;
  n = 10;
  two = 2;
  alpha = nu;
  fnu = nu;

  p = c_prod(nref, k0);
  /* p >= beta */
  state = c_comparege(p, beta);
  if (state)
    {
      p1 = c_prod(nref, k0);
      pow_ci(&p2, &p1, &two);
      pow_ci(&p, &beta, &two);
      p = c_sub(p2, p);
      c_sqrt(&kappa, &p);

      p1 = cscalar_prod(beta, nu);
      pow_ci(&p2, &kappa, &two);
      p2 = cscalar_prod(p2, r1);
      c_div(&p, &p1, &p2);
      alpha1 = p;

      p1 = cscalar_prod(beta, nu);
      pow_ci(&p2, &kappa, &two);
      p2 = cscalar_prod(p2, r2);
      c_div(&p, &p1, &p2);      
      alpha2 = p;

      c_div(&rho, &k0, &kappa);

      pow_ci(&p1, &nref, &two);      
      p2 = c_prod(k0, p1);
      c_div(&p, &p2, &kappa);
      sigma = p;

      z1 = cscalar_prod(kappa, r1);
      z2 = cscalar_prod(kappa, r2);

      cbesj(&z2, &alpha, &kode, &n, del1bjy, &nz, &ierr);
      cbesyp(&z1, &fnu, &kode, &n, del1bypy, &nz, cwrk, &ierr);
      cbesjp(&z1, &alpha, &kode, &n, del1bjpy, &nz, &ierr);
      cbesy(&z2, &fnu, &kode, &n, del1byy, &nz, cwrk, &ierr);

      p1 = c_prod(del1bjy[0], del1bypy[0]);
      p2 = c_prod(del1bjpy[0], del1byy[0]);
      p = c_sub(p1, p2);
      del1 = p;
      
      cbesj(&z1, &alpha, &kode, &n, del2bjy, &nz, &ierr);
      cbesyp(&z2, &fnu, &kode, &n, del2bypy, &nz, cwrk, &ierr);
      cbesjp(&z2, &alpha, &kode, &n, del2bjpy, &nz, &ierr);
      cbesy(&z1, &fnu, &kode, &n, del2byy, &nz, cwrk, &ierr);

      p1 = c_prod(del2bjy[0], del2bypy[0]);
      p2 = c_prod(del2bjpy[0], del2byy[0]);
      p = c_sub(p1, p2);
      del2 = p;

      cbesj(&z2, &alpha, &kode, &n, del3bjy2, &nz, &ierr);
      cbesy(&z1, &fnu, &kode, &n, del3byy1, &nz, cwrk, &ierr);
      cbesj(&z1, &alpha, &kode, &n, del3bjy1, &nz, &ierr);
      cbesy(&z2, &fnu, &kode, &n, del3byy2, &nz, cwrk, &ierr);

      p1 = c_prod(del3bjy2[0], del3byy1[0]);
      p2 = c_prod(del3bjy1[0], del3byy2[0]);
      p = c_sub(p1, p2);
      del3 = p;

      cbesjp(&z2, &alpha, &kode, &n, del4bjpy2, &nz, &ierr);
      cbesyp(&z1, &fnu, &kode, &n, del4bypy1, &nz, cwrk, &ierr);
      cbesjp(&z1, &alpha, &kode, &n, del4bjpy1, &nz, &ierr);
      cbesyp(&z2, &fnu, &kode, &n, del4bypy2, &nz, cwrk, &ierr);

      p1 = c_prod(del4bjpy2[0], del4bypy1[0]);
      p2 = c_prod(del4bjpy1[0], del4bypy2[0]);
      p = c_sub(p1, p2);
      del4 = p;

      lay[0][0] = del1;

      p1 = c_prod(alpha1, del3);
      c_div(&p, &p1, &sigma);
      lay[0][1] = p;
      lay[0][2] = czero;

      c_div(&p1, &del3, &sigma);
      p = cscalar_prod(p1, (-1.0));
      lay[0][3] = p;
      
      p1 = c_prod(alpha1, del3);
      c_div(&p, &p1, &rho);
      lay[1][0] = p;
      lay[1][1] = del1;

      c_div(&p1, &del3, &rho);
      p = cscalar_prod(p1, (-1.0));
      lay[1][2] = p;
      lay[1][3] = czero;

      p1 = c_prod(alpha2, del1);
      p2 = c_prod(alpha1, del2);
      p = c_sub(p1, p2);
      lay[2][0] = p;

      p1 = c_prod(alpha1, alpha2);
      p2 = c_prod(p1, del3);
      c_div(&p, &p2, &sigma);
      p1 = c_prod(rho, del4);
      p = c_add(p, p1);
      lay[2][1] = p;
      lay[2][2] = del2;

      p1 = c_prod(alpha2, del3);
      c_div(&p2, &p1, &sigma);
      p = cscalar_prod(p2, (-1.0));
      lay[2][3] = p;
      
      p1 = c_prod(alpha1, alpha2);
      p1 = c_prod(p1, del3);
      c_div(&p2, &p1, &rho);
      p1 = c_prod(sigma, del4);
      p = c_add(p2, p1);
      lay[3][0] = p;

      p1 = c_prod(alpha2, del1);
      p2 = c_prod(alpha1, del2);
      p = c_sub(p1, p2);
      lay[3][1] = p;

      p1 = c_prod(alpha2, del3);
      c_div(&p2, &p1, &rho);
      p = cscalar_prod(p2, (-1.0));
      lay[3][2] = p;

      lay[3][3] = del2;

      for (i=0;i<4;i++)
	{
	  for (j=0;j<4;j++)
	    {
	      p1 = cscalar_prod(z1, PI);
	      p = cscalar_prod(p1, (0.5));
	      ftm[i][j] = c_prod(lay[i][j],p);
	    }
	}
      return (1);
    }
  
  else
    {
      pow_ci(&p1, &beta, &two);
      p2 = c_prod(nref, k0);
      pow_ci(&p, &p2, &two);
      p = c_sub(p1, p);
      c_sqrt(&gamma, &p);

      p1 = cscalar_prod(beta, nu);
      pow_ci(&p2, &gamma, &two);
      p2 = cscalar_prod(p2, r1);
      c_div(&p, &p1, &p2);
      alpha1 = p;

      p1 = cscalar_prod(beta, nu);
      pow_ci(&p2, &gamma, &two);
      p2 = cscalar_prod(p2, r2);
      c_div(&p, &p1, &p2);      
      alpha2 = p;

      c_div(&rho, &k0, &gamma);

      pow_ci(&p1, &nref, &two);
      p1 = c_prod(p1, k0);
      c_div(&sigma, &p1, &gamma);

      z1 = cscalar_prod(gamma, r1);
      z2 = cscalar_prod(gamma, r2);

      cbesip(&z1, &alpha, &kode, &n, del1bipy, &nz, &ierr);
      cbesk(&z2, &fnu, &kode, &n, del1bky, &nz, &ierr);
      cbesi(&z2, &alpha, &kode, &n, del1biy, &nz, &ierr);
      cbeskp(&z1, &fnu, &kode, &n, del1bkpy, &nz, &ierr);

      p1 = c_prod(del1bipy[0], del1bky[0]);
      p2 = c_prod(del1biy[0], del1bkpy[0]);
      del1 = c_sub(p1, p2);
      
      cbesip(&z2, &alpha, &kode, &n, del2bipy, &nz, &ierr);
      cbesk(&z1, &fnu, &kode, &n, del2bky, &nz, &ierr);
      cbesi(&z1, &alpha, &kode, &n, del2biy, &nz, &ierr);
      cbeskp(&z2, &fnu, &kode, &n, del2bkpy, &nz, &ierr);

      p1 = c_prod(del2bipy[0], del2bky[0]);
      p2 = c_prod(del2biy[0], del2bkpy[0]);
      del2 = c_sub(p1, p2);

      cbesi(&z1, &alpha, &kode, &n, del3biy2, &nz, &ierr);
      cbesk(&z2, &fnu, &kode, &n, del3bky1, &nz, &ierr);
      cbesi(&z2, &alpha, &kode, &n, del3biy1, &nz, &ierr);
      cbesk(&z1, &fnu, &kode, &n, del3bky2, &nz, &ierr);

      p1 = c_prod(del3biy2[0], del3bky1[0]);
      p2 = c_prod(del3biy1[0], del3bky2[0]);
      del3 = c_sub(p1, p2);

      cbesip(&z1, &alpha, &kode, &n, del4bipy2, &nz, &ierr);
      cbeskp(&z2, &fnu, &kode, &n, del4bkpy1, &nz, &ierr);
      cbesip(&z2, &alpha, &kode, &n, del4bipy1, &nz, &ierr);
      cbeskp(&z1, &fnu, &kode, &n, del4bkpy2, &nz, &ierr);

      p1 = c_prod(del4bipy2[0], del4bkpy1[0]);
      p2 = c_prod(del4bipy1[0], del4bkpy2[0]);
      del4 = c_sub(p1, p2);

      lay[0][0] = del1;
      
      p1 = c_prod(alpha1, del3);
      c_div(&p, &p1, &sigma);
      lay[0][1] = p;
      lay[0][2] = czero;

      c_div(&p, &del3, &sigma);
      lay[0][3] = p;
      
      p1 = c_prod(alpha1, del3);
      c_div(&p, &p1, &rho);
      lay[1][0] = p;
      lay[1][1] = del1;

      c_div(&p, &del3, &rho);
      lay[1][2] = p;
      lay[1][3] = czero;

      p1 = c_prod(alpha2, del1);
      p2 = c_prod(alpha1, del2);
      p = c_sub(p2, p1);
      lay[2][0] = p;

      p1 = c_prod(alpha1, alpha2);
      p1 = c_prod(p1, del3);
      c_div(&p2, &p1, &sigma);
      p2 = cscalar_prod(p2, (-1.0));
      p1 = c_prod(rho, del4);
      p = c_sub(p2, p1);
      lay[2][1] = p;
      lay[2][2] = del2;

      p1 = c_prod(alpha2, del3);
      c_div(&p2, &p1, &sigma);
      p = cscalar_prod(p2, (-1.0));
      lay[2][3] = p;
      
      p1 = c_prod(alpha1, alpha2);
      p1 = c_prod(p1, del3);
      c_div(&p2, &p1, &rho);
      p2 = cscalar_prod(p2, (-1.0));
      p1 = c_prod(sigma, del4);
      p = c_sub(p2, p1);
      lay[3][0] = p;

      p1 = c_prod(alpha2, del1);
      p2 = c_prod(alpha1, del2);
      p = c_sub(p2, p1);
      lay[3][1] = p;

      p1 = c_prod(alpha2, del3);
      c_div(&p2, &p1, &rho);
      p = cscalar_prod(p2, (-1.0));
      lay[3][2] = p;
      lay[3][3] = del2;

      for (i=0;i<4;i++)
	{
	  for (j=0;j<4;j++)
	    {
	      p = c_prod(lay[i][j], z1);
	      ftm[i][j] = p;
	    }
	}
      return (-1);
    }
}
