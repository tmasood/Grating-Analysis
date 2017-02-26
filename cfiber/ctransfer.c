#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"
 
#define PI M_PI

#define SIZE 4

extern int cbesyp(complex*, real*, integer*, integer*, complex*,
		  integer*, complex*, integer*);
extern int cbesy(complex*, real*, integer*, integer*, complex*,
		 integer*, complex*, integer*);
extern int cbesjp(complex*, real*, integer*, integer*, complex*,
		  integer*, integer*);
extern int cbesj(complex*, real*, integer*, integer*, complex*,
		 integer*, integer*);
extern int cbesip(complex*, real*, integer*, integer*, complex*,
		  integer*, integer*);
extern int cbesk(complex*, real*, integer*, integer*, complex*,
		 integer*, integer*);
extern int cbesi(complex*, real*, integer*, integer*, complex*,
		 integer*, integer*);
extern int cbeskp(complex*, real*, integer*, integer*, complex*,
		  integer*, integer*);
extern int cmatmul(complex*, complex*, complex*);
 
/* find Coefficient transfer matrix for a given layer */
int ctransfer(complex beta, real nu, complex n2, complex n1,
	      real r1, complex k0, complex MatP[][SIZE])
{
  complex M_product[SIZE][SIZE];
  complex Mat2[SIZE][SIZE];
  complex Mat1[SIZE][SIZE];
  complex kappa2;
  complex alpha2;
  complex rho2;
  complex z2yp[10], z2y[10], z2jp[10], z2j[10];
  complex z2i[10], z2ip[10], z2k[10], z2kp[10], z1j[10];
  complex z1jp[10], z1yp[10], z1y[10];
  complex z1k[10], z1kp[10], z1ip[10], z1i[10];
  complex cwrk[10];
  complex kappa1;
  complex alpha1;
  real alpha;
  complex rho1;
  complex sigma1;
  complex sigma2;
  complex gamma2;
  complex gamma1;
  complex z2, z1;
  complex p1, p2, p;
  complex czero;
  integer kode, ierr;
  integer n, nz;
  integer two;
  int i,j;
  int state;
  n = 10;

  czero.r = 0.0;
  czero.i = 0.0;

  two = 2;
  kappa1 = czero;
  kode = 1;
  alpha = nu;

  p = c_prod(n2, k0);
  state = c_comparege(p, beta);
  /* p >= beta */
  if (state)
    {
      p1 = c_prod(n2, k0);
      pow_ci(&p2, &p1, &two);
      pow_ci(&p, &beta, &two);
      p = c_sub(p2, p);
      c_sqrt(&kappa2, &p);      

      p1 = cscalar_prod(beta, nu);
      pow_ci(&p2, &kappa2, &two);
      p2 = cscalar_prod(p2, r1);
      c_div(&alpha2, &p1, &p2);

      c_div(&rho2, &k0, &kappa2);

      pow_ci(&p1, &n2, &two);      
      p2 = c_prod(k0, p1);
      c_div(&sigma2, &p2, &kappa2);

      z2 = cscalar_prod(kappa2, r1);

      cbesyp(&z2, &nu, &kode, &n, z2yp, &nz, cwrk, &ierr);
      cbesy(&z2, &nu, &kode, &n, z2y, &nz, cwrk, &ierr);
      cbesjp(&z2, &alpha, &kode, &n, z2jp, &nz, &ierr);
      cbesj(&z2, &alpha, &kode, &n, z2j, &nz, &ierr);

      Mat2[0][0] = z2yp[0];

      p1 = c_prod(alpha2, z2y[0]);
      c_div(&p, &p1, &sigma2);
      Mat2[0][1] = p;
      Mat2[0][2] = czero;

      c_div(&p, &z2y[0], &sigma2);
      p = cscalar_prod(p, (-1.0));
      Mat2[0][3] = p;

      p = cscalar_prod(z2jp[0], (-1.0));
      Mat2[1][0] = p;

      p1 = c_prod(alpha2, z2j[0]);
      c_div(&p, &p1, &sigma2);
      p = cscalar_prod(p, (-1.0));
      Mat2[1][1] = p;
      Mat2[1][2] = czero;

      c_div(&p, &z2j[0], &sigma2);
      Mat2[1][3] = p;

      p1 = c_prod(alpha2, z2y[0]);
      c_div(&p, &p1, &rho2);
      Mat2[2][0] = p;
      Mat2[2][1] = z2yp[0];

      c_div(&p1, &z2y[0], &rho2);
      p = cscalar_prod(p1, (-1.0));
      Mat2[2][2] = p;
      Mat2[2][3] = czero;

      p1 = c_prod(alpha2, z2j[0]);
      c_div(&p, &p1, &rho2);
      p = cscalar_prod(p, (-1.0));      
      Mat2[3][0] = p;

      p = cscalar_prod(z2jp[0], (-1.0));
      Mat2[3][1] = p;

      c_div(&p, &z2j[0], &rho2);
      Mat2[3][2] = p;
      Mat2[3][3] = czero;

      for (i=0; i<4; i++)
	{
	  for (j=0; j<4; j++)
	    {
	      p1 = cscalar_prod(z2, (PI/2));
	      p = c_prod(p1, Mat2[i][j]);
	      Mat2[i][j] = p;
	    }
	}
    }
  else
    {
      p1 = c_prod(n2, k0);
      pow_ci(&p2, &p1, &two);
      pow_ci(&p, &beta, &two);
      p = c_sub(p, p2);
      c_sqrt(&gamma2, &p);      

      p1 = cscalar_prod(beta, nu);
      pow_ci(&p2, &gamma2, &two);
      p2 = cscalar_prod(p2, r1);
      c_div(&alpha2, &p1, &p2);

      c_div(&rho2, &k0, &gamma2);

      pow_ci(&p1, &n2, &two);      
      p2 = c_prod(k0, p1);
      c_div(&sigma2, &p2, &gamma2);

      z2 = cscalar_prod(gamma2, r1);

      cbesi(&z2, &alpha, &kode, &n, z2i, &nz, &ierr);
      cbesip(&z2, &alpha, &kode, &n, z2ip, &nz, &ierr);
      cbesk(&z2, &nu, &kode, &n, z2k, &nz, &ierr);
      cbeskp(&z2, &nu, &kode, &n, z2kp, &nz, &ierr);

      Mat2[0][0] = z2ip[0];

      p1 = c_prod(alpha2, z2i[0]);
      c_div(&p, &p1, &sigma2);
      Mat2[0][1] = p;
      Mat2[0][2] = czero;

      c_div(&p, &z2i[0], &sigma2);
      Mat2[0][3] = p;

      p = cscalar_prod(z2kp[0], (-1.0));
      Mat2[1][0] = p;

      p1 = c_prod(alpha2, z2k[0]);
      c_div(&p2, &p1, &sigma2);
      p = cscalar_prod(p2, (-1.0));
      Mat2[1][1] = p;
      Mat2[1][2] = czero;

      c_div(&p1, &z2k[0], &sigma2);
      p = cscalar_prod(p1, (-1.0));
      Mat2[1][3] = p;

      p1 = c_prod(alpha2, z2i[0]);
      c_div(&p, &p1, &rho2);
      Mat2[2][0] = p;
      Mat2[2][1] = z2ip[0];

      c_div(&p, &z2i[0], &rho2);
      Mat2[2][2] = p;
      Mat2[2][3] = czero;

      p1 = c_prod(alpha2, z2k[0]);
      c_div(&p, &p1, &rho2);
      p = cscalar_prod(p, (-1.0));
      Mat2[3][0] = p;

      p = cscalar_prod(z2kp[0], (-1.0));
      Mat2[3][1] = p;

      c_div(&p, &z2k[0], &rho2);
      p = cscalar_prod(p, (-1.0));
      Mat2[3][2] = p;
      Mat2[3][3] = czero;
      
      for (i=0; i<4; i++)
	{
	  for (j=0; j<4; j++)
	    {
	      Mat2[i][j] = c_prod(z2, Mat2[i][j]);
	    }
	}
    }

  p = c_prod(n1, k0);
  state = c_comparege(p, beta);
  /* p >= beta */
  if (state)
    {
      p1 = c_prod(n1, k0);
      pow_ci(&p2, &p1, &two);
      pow_ci(&p, &beta, &two);
      p = c_sub(p2, p);
      c_sqrt(&kappa1, &p);
      
      p1 = cscalar_prod(beta, nu);
      pow_ci(&p2, &kappa1, &two);
      p2 = cscalar_prod(p2, r1);
      c_div(&alpha1, &p1, &p2);

      c_div(&rho1, &k0, &kappa1);

      pow_ci(&p1, &n1, &two);      
      p2 = c_prod(k0, p1);
      c_div(&sigma1, &p2, &kappa1);

      z1 = cscalar_prod(kappa1, r1);

      cbesj(&z1, &alpha, &kode, &n, z1j, &nz, &ierr);
      cbesjp(&z1, &alpha, &kode, &n, z1jp, &nz, &ierr);
      cbesyp(&z1, &nu, &kode, &n, z1yp, &nz, cwrk, &ierr);
      cbesy(&z1, &nu, &kode, &n, z1y, &nz, cwrk, &ierr);

      Mat1[0][0] = z1j[0];
      Mat1[0][1] = z1y[0];
      Mat1[0][2] = czero;
      Mat1[0][3] = czero;
      
      Mat1[1][0] = czero;
      Mat1[1][1] = czero;
      Mat1[1][2] = z1j[0];
      Mat1[1][3] = z1y[0];

      Mat1[2][0] = c_prod(alpha1,z1j[0]);
      Mat1[2][1] = c_prod(alpha1,z1y[0]);
      Mat1[2][2] = c_prod(rho1, z1jp[0]);
      Mat1[2][3] = c_prod(rho1, z1yp[0]);

      Mat1[3][0] = c_prod(sigma1, z1jp[0]);
      Mat1[3][1] = c_prod(sigma1, z1yp[0]);
      Mat1[3][2] = c_prod(alpha1, z1j[0]);
      Mat1[3][3] = c_prod(alpha1, z1y[0]);
    }
  else
    {
      p1 = c_prod(n1, k0);
      pow_ci(&p2, &p1, &two);
      pow_ci(&p, &beta, &two);
      p = c_sub(p, p2);
      c_sqrt(&gamma1, &p);    

      p1 = cscalar_prod(beta, nu);
      pow_ci(&p2, &gamma1, &two);
      p2 = cscalar_prod(p2, r1);
      c_div(&alpha1, &p1, &p2);

      c_div(&rho1, &k0, &gamma1);

      pow_ci(&p1, &n1, &two);      
      p2 = c_prod(k0, p1);
      c_div(&sigma1, &p2, &gamma1);

      z1 = cscalar_prod(gamma1, r1);

      cbesi(&z1, &alpha, &kode, &n, z1i, &nz, &ierr);
      cbesip(&z1, &alpha, &kode, &n, z1ip, &nz, &ierr);
      cbesk(&z1, &nu, &kode, &n, z1k, &nz, &ierr);
      cbeskp(&z1, &nu, &kode, &n, z1kp, &nz, &ierr);

      Mat1[0][0] = z1k[0];
      Mat1[0][1] = z1i[0];
      Mat1[0][2] = czero;
      Mat1[0][3] = czero;
      
      Mat1[1][0] = czero;
      Mat1[1][1] = czero;
      Mat1[1][2] = z1k[0];
      Mat1[1][3] = z1i[0];

      p1 = c_prod(alpha1, z1k[0]);
      p = cscalar_prod(p1, (-1.0));
      Mat1[2][0] = p;

      p1 = c_prod(alpha1, z1i[0]);
      p = cscalar_prod(p1, (-1.0));
      Mat1[2][1] = p;

      p1 = c_prod(rho1, z1kp[0]);
      p = cscalar_prod(p1, (-1.0));
      Mat1[2][2] = p;

      p1 = c_prod(rho1, z1ip[0]);
      p = cscalar_prod(p1, (-1.0));
      Mat1[2][3] = p;

      p1 = c_prod(sigma1, z1kp[0]);
      p = cscalar_prod(p1, (-1.0));
      Mat1[3][0] = p;

      p1 = c_prod(sigma1, z1ip[0]);
      p = cscalar_prod(p1, (-1.0));
      Mat1[3][1] = p;

      p1 = c_prod(alpha1, z1k[0]);
      p = cscalar_prod(p1, (-1.0));
      Mat1[3][2] = p;

      p1 = c_prod(alpha1, z1i[0]);
      p = cscalar_prod(p1, (-1.0));
      Mat1[3][3] = p;
    }

  cmatmul(Mat2, Mat1, M_product);
  for (i=0; i<4; i++)
    {
      for (j=0; j<4; j++)
	{
	  MatP[i][j] = M_product[i][j];
	}
    }

  p = c_prod(n2, k0);
  /* p >= beta */
  state = c_comparege(p, beta);

  if (state)
    {
      return (1);
    }
  else
    {
      return (-1);
    }
}
