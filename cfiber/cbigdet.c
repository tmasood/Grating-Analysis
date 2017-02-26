#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"

#define PI M_PI
#define SIZE 10

extern int cproduct(int, complex, real, complex, complex*, real*, complex*);
extern int cbesj(complex*, real*, integer*, integer*, complex*,
		 integer*, integer*);
extern int cbesjp(complex*, real*, integer*, integer*, complex*, integer*,
		  integer*);
extern int cbesi(complex*, real*, integer*, integer*, complex*, integer*,
		 integer*);
extern int cbesip(complex*, real*, integer*, integer*, complex*, integer*,
		  integer*);
extern int cmatmul(complex*, complex*, complex*);
extern int cbesk(complex*, real*, integer*, integer*, complex*, integer*,
		 integer*);
extern int cbeskp(complex*, real*, integer*, integer*, complex*, integer*,
		  integer*);
extern int ccaldet(complex*, complex*);

/* find total matrix product as function of beta */
int cbigdet(int nbound, complex *wbeta, real nu, complex k0,
	   complex *index, real *wrb, complex *bigdetc)
{
  complex M_product[4][4];
  complex temp[4][4];
  complex CharMat[4][4];
  complex M0[4][4];
  
  complex kappa;
  complex gamma;
  complex gamma2;
  complex p1, p2, p;
  complex foo, bar;
  complex z1;
  complex z2;
  integer n = 10;
  integer nz;
  integer kode;
  integer ierr;
  complex cbesjm0[SIZE], cbesjpm0[SIZE], cbesim0[SIZE], cbesipm0[SIZE];
  complex cbeskcm[SIZE], cbeskpcm[SIZE];
  complex temp1;
  complex czero;
  complex cb;
  real vflag;

  integer two;

  int m;
  int i, j;
  int state;

  kode = 1;
  m = 0;
  two = 2;

  czero.r = 0.0;
  czero.i = 0.0;

  while (m < 100)
    {
      cproduct(nbound, wbeta[m], nu, k0, index, wrb, M_product);
      for (i=0; i<4; i++)
	{
	  for (j=0; j<4; j++)
	    {
	      temp[i][j] = M_product[i][j];
	    }
	}
      pow_ci(&p1,&wbeta[m],&two);
      pow_ci(&p1,&wbeta[m],&two);
      cb = c_prod(index[nbound],k0);
      pow_ci(&p2,&cb,&two);
      p = c_sub(p1, p2);
      c_sqrt(&gamma2, &p);
      z2 = cscalar_prod(gamma2, wrb[nbound]);
      cb = c_prod(index[0],k0);
      pow_ci(&p2,&cb,&two);
      pow_ci(&p1,&wbeta[m],&two);
      p = c_sub(p2, p1);
      state = c_comparege(p, czero);

      if (state)
	{
	  vflag = 1.0;
	  c_sqrt(&kappa, &p);
	  z1 = cscalar_prod(kappa, wrb[1]);
	  cbesj(&z1, &nu, &kode, &n, cbesjm0, &nz, &ierr);
	  cbesjp(&z1, &nu, &kode, &n, cbesjpm0, &nz, &ierr);

	  M0[0][0] = cbesjm0[0];
	  M0[0][1] = czero;
	  M0[0][2] = czero;
	  M0[0][3] = czero;
	  
	  M0[1][0] = czero;
	  M0[1][1] = cbesjm0[0];
	  M0[1][2] = czero;
	  M0[1][3] = czero;
	  
	  p1 = c_prod(kappa,z1);
	  p2 = cscalar_prod(wbeta[m], nu);
	  c_div(&p, &p2, &p1);
	  p = c_prod(p, cbesjm0[0]);

	  M0[2][0] = p;
	  c_div(&p1, &k0, &kappa);
	  p = c_prod(p1, cbesjpm0[0]);
	  M0[2][1] = p;
	  M0[2][2] = czero;
	  M0[2][3] = czero;

	  pow_ci(&p1,&index[0],&two);
	  p2 = c_prod(p1, k0);
	  c_div(&p, &p2, &kappa);
	  p = c_prod(p, cbesjpm0[0]);
	  M0[3][0] = p;

	  p1 = cscalar_prod(wbeta[m],nu);
	  p2 = c_prod(kappa, z1);
	  c_div(&p, &p1, &p2);
	  p = c_prod(p, cbesjm0[0]);
	  M0[3][1] = p;
	  M0[3][2] = czero;
	  M0[3][3] = czero;
	}
      else
	{
	  vflag = -1.0;
	  pow_ci(&p1,&wbeta[m],&two);
	  cb = c_prod(index[0],k0);
	  pow_ci(&p2,&cb,&two);
	  p.r = p1.r - p2.r;
	  p.i = p1.i - p2.i;	  
	  c_sqrt(&gamma, &p);
	  z1 = cscalar_prod(gamma, wrb[1]);
	  cbesi(&z1, &nu, &kode, &n, cbesim0, &nz, &ierr);
	  cbesip(&z1, &nu, &kode, &n, cbesipm0, &nz, &ierr);
	  M0[0][0] = cbesim0[0];
	  M0[0][1] = czero;
	  M0[0][2] = czero;
	  M0[0][3] = czero;
	  
	  M0[1][0] = czero;
	  M0[1][1] = cbesim0[0];
	  M0[1][2] = czero;
	  M0[1][3] = czero;

	  p1 = c_prod(gamma,z1);
	  p2 = cscalar_prod(wbeta[m], nu);
	  c_div(&p, &p2, &p1);
	  p = c_prod(p, cbesim0[0]);
	  p = cscalar_prod(p, (-1.0));

	  M0[2][0] = p;
	  c_div(&p1, &k0, &gamma);
	  p2 = c_prod(p1, cbesim0[0]);
	  p = cscalar_prod(p2, (-1.0));
	  M0[2][1] = p;
	  M0[2][2] = czero;
	  M0[2][3] = czero;
	  
	  pow_ci(&p1, &index[0], &two);
	  p2 = c_prod(p1, k0);
	  c_div(&p, &p2, &gamma);
	  p = c_prod(p, cbesipm0[0]);
	  p = cscalar_prod(p, (-1.0));
	  M0[3][0] = p;

	  p1 = cscalar_prod(wbeta[m], nu);
	  p2 = c_prod(gamma, z1);
	  c_div(&p, &p1, &p2);
	  p = c_prod(p, cbesim0[0]);
	  p = cscalar_prod(p, (-1.0));
	  M0[3][1] = p;
	  M0[3][2] = czero;
	  M0[3][3] = czero;
	}
      cmatmul(temp, M0, M_product);
      cbesk(&z2, &nu, &kode, &n, cbeskcm, &nz, &ierr);
      cbeskp(&z2, &nu, &kode, &n, cbeskpcm, &nz, &ierr);
      CharMat[0][0] = M_product[0][0];
      CharMat[0][1] = M_product[0][1];
      p = cscalar_prod(cbeskcm[0],(-1.0));
      CharMat[0][2] = p;
      CharMat[0][3] = czero;

      CharMat[1][0] = M_product[1][0];
      CharMat[1][1] = M_product[1][1];
      CharMat[1][2] = czero;
      p = cscalar_prod(cbeskcm[0], (-1.0));
      CharMat[1][3] = p;
      
      CharMat[2][0] = M_product[2][0];
      CharMat[2][1] = M_product[2][1];
      p1 = cscalar_prod(wbeta[m], nu);
      p2 = c_prod(gamma2, z2);
      c_div(&p, &p1, &p2);
      p = c_prod(p, cbeskcm[0]);
      CharMat[2][2] = p;
      c_div(&p1, &k0, &gamma2);
      p = c_prod(p1, cbeskpcm[0]);
      CharMat[2][3] = p;

      CharMat[3][0] = M_product[3][0];
      CharMat[3][1] = M_product[3][1];
      
      pow_ci(&p1,&index[nbound],&two);
      p2 = c_prod(p1, k0);
      c_div(&p, &p2, &gamma2);
      p = c_prod(p, cbeskpcm[0]);
      CharMat[3][2] = p;

      p1 = cscalar_prod(wbeta[m], nu);
      p2 = c_prod(gamma2, z2);
      c_div(&p, &p1, &p2);
      p = c_prod(p, cbeskcm[0]);
      CharMat[3][3] = p;

      if (nu == 0.0)
	{
	  vflag = 1.0;
	}
      ccaldet(CharMat, &temp1); 
      bigdetc[m] = temp1;
      m++;
    }
  return (0);  
}
