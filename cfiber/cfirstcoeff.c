#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"
 
#define PI M_PI
#define SIZE 4

extern int cproduct(int, complex, real, complex, complex*, real*, complex*);
extern int cbesj(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbesjp(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbesip(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbesk(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbesi(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbeskp(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cmatmul(complex*, complex*, complex*);
extern int cmatsolv(complex*, complex*, complex*);
/* find first set of eigenmode coefficients */
int cfirstcoeff(int nbound, complex betaR, real nu, complex k0, complex *index,
	       real *wrb, complex *M_x)
{
  complex M_product[SIZE][SIZE];
  complex temp[SIZE][SIZE];
  complex EigmodeMat[SIZE - 1][SIZE - 1];
  complex EigVect[SIZE - 1];
  complex M0[SIZE][SIZE];
  complex besjm0[10], besjpm0[10],besim0[10], besipm0[10];
  complex beskcm[10], beskpcm[10];
  complex gamma2;
  complex gamma;
  complex z2, z1;
  complex kappa;
  complex p1, p2, p;
  complex czero;
 
  int i,j;
  int state;
  integer n = 10;
  integer nz;
  integer  kode;
  integer ierr;
  integer two;

  two = 2;
  kode = 1;
  
  cproduct(nbound, betaR, nu, k0, index, wrb, M_product);
  for (i=0; i<SIZE; i++)
    {
      for (j=0; j<SIZE; j++)
	{
	  temp[i][j] = M_product[i][j];
	}
    }

  p1 = c_prod(index[nbound], k0);
  pow_ci(&p2, &p1, &two);
  pow_ci(&p, &betaR, &two);
  p = c_sub(p, p2);
  c_sqrt(&gamma2, &p);    

  z2 = cscalar_prod(gamma2, wrb[nbound]);

  pow_ci(&p1, &betaR, &two);
  p2 = c_prod(index[0], k0);
  pow_ci(&p, &p2, &two);
  p = c_sub(p, p1);
  /* p >= 0 */
  state = c_comparege(p, czero);

  if (state)
    {
      p1 = c_prod(index[0], k0);
      pow_ci(&p2, &p1, &two);
      pow_ci(&p1, &betaR, &two);
      p = c_sub(p2, p1);
      c_sqrt(&kappa, &p);

      z1 = cscalar_prod(kappa, wrb[1]);

      cbesj(&z1, &nu, &kode, &n, besjm0, &nz, &ierr);
      cbesjp(&z1, &nu, &kode, &n, besjpm0, &nz, &ierr);
      M0[0][0] = besjm0[0];
      M0[0][1] = czero;
      M0[0][2] = czero;
      M0[0][3] = czero;
      
      M0[1][0] = czero;
      M0[1][1] = besjm0[0];
      M0[1][2] = czero;
      M0[1][3] = czero;
      
      p1 = cscalar_prod(betaR, nu);
      p2 = c_prod(kappa, z1);
      c_div(&p, &p1, &p2);
      p = c_prod(p, besjm0[0]);
      M0[2][0] = p;

      c_div(&p1, &k0, &kappa);
      p = c_prod(p1, besjpm0[0]);
      M0[2][1] = p;
      M0[2][2] = czero;
      M0[2][3] = czero;
 
      pow_ci(&p1, &index[0], &two);
      p2 = c_prod(p1, k0);
      c_div(&p1, &p2, &kappa);
      p = c_prod(p1, besjpm0[0]);
      M0[3][0] = p;

      p1 = cscalar_prod(betaR, nu);
      p2 = c_prod(kappa, z1);
      c_div(&p, &p1, &p2);
      p = c_prod(p, besjm0[0]);
      M0[3][1] = p;
      M0[3][2] = czero;
      M0[3][3] = czero;
    }
  else
    {
      pow_ci(&p1, &betaR, &two);
      p2 = c_prod(index[0], k0);
      pow_ci(&p, &p2, &two);
      p = c_sub(p1, p);
      c_sqrt(&gamma, &p);

      z1 = cscalar_prod(gamma, wrb[1]);

      cbesi(&z1, &nu, &kode, &n, besim0, &nz, &ierr);
      cbesip(&z1, &nu, &kode, &n, besipm0, &nz, &ierr);
      M0[0][0] = besim0[0];
      M0[0][1] = czero;
      M0[0][2] = czero;
      M0[0][3] = czero;
          
      M0[1][0] = czero;
      M0[1][1] = besim0[0];
      M0[1][2] = czero;
      M0[1][3] = czero;
      
      p1 = cscalar_prod(betaR, nu);
      p2 = c_prod(gamma, z1);
      c_div(&p, &p1, &p2);
      p = c_prod(p, besim0[0]);
      p = cscalar_prod(p, (-1.0));
      M0[2][0] = p;

      c_div(&p1, &k0, &gamma);
      p2 = c_prod(p1, besipm0[0]);
      p = cscalar_prod(p2, (-1.0));
      M0[2][1] = p;
      M0[2][2] = czero;
      M0[2][3] = czero;
          
      pow_ci(&p1, &index[0], &two);
      p2 = c_prod(p1, k0);
      c_div(&p1, &p2, &gamma);
      p = c_prod(p1, besipm0[0]);
      p = cscalar_prod(p, (-1.0));
      M0[3][0] = p;

      p1 = cscalar_prod(betaR, nu);
      p2 = c_prod(gamma, z1);
      c_div(&p, &p1, &p2);
      p = c_prod(p, besim0[0]);
      p = cscalar_prod(p, (-1.0));
      M0[3][1] = p;
      M0[3][2] = czero;
      M0[3][3] = czero;
    }

  cmatmul(temp, M0, M_product);
  cbesk(&z2, &nu, &kode, &n, beskcm, &nz, &ierr);
  cbeskp(&z2, &nu, &kode, &n, beskpcm, &nz, &ierr);

      EigmodeMat[0][0] = M_product[0][1];
      p = cscalar_prod(beskcm[0], (-1.0));
      EigmodeMat[0][1] = p;
      EigmodeMat[0][2] = czero;
 
      EigmodeMat[1][0] = M_product[1][1];
      EigmodeMat[1][1] = czero;
      
      p = cscalar_prod(beskcm[0], (-1.0));
      EigmodeMat[1][2] = p;
      
      EigmodeMat[2][0] = M_product[2][1];
      
      p1 = cscalar_prod(betaR, nu);
      p2 = c_prod(gamma2, z2);
      c_div(&p, &p1, &p2);
      p = c_prod(p, beskcm[0]);      
      EigmodeMat[2][1] = p;

      c_div(&p1, &k0, &gamma2);
      p = c_prod(p1, beskpcm[0]);
      EigmodeMat[2][2] = p;

      /*     for (i=0;i<3;i++)
	{
	  for (j=0;j<3;j++)
	    {
	      printf("EigenModeMat[%d][%d] = %f \n",i,j,EigmodeMat[i][j]);
	    }
	}
	*/
 
      p = cscalar_prod(M_product[0][0], (-1.0));
      EigVect[0] = p;

      p = cscalar_prod(M_product[1][0], (-1.0));
      EigVect[1] = p;

      p = cscalar_prod(M_product[2][0], (-1.0));
      EigVect[2] = p;

      /*	  for (j=0;j<3;j++)
	    {
	      printf("EigenVect[%d] = %f \n",j,EigVect[j]);
	    }
	    */
      cmatsolv(EigmodeMat, EigVect, M_x);
      return (0);
}
