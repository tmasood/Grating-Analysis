#include <stdio.h>
#include <math.h>
#include "f2c.h"

extern int besj_();
extern int cbesj_();
extern int besjp();
extern int besy_();
extern int besi_();

int main()
{
  complex z__;
  real fnu;
  integer kode, n;
  complex cy;
  integer nz, ierr; 
  real x, alpha;
  real y[10];
  real oldx;
  real oldzr;

  real beta, k0;
  real index;
  real wrb0, wrb1;
  real Aindex[10];
  real Awrb[10];
  real wbeta[100];
  real ftm[4][4];
  real M_product[4][4];
  int i,j;
  int nbound;

  real bigdetc[100];

  /*  printf("\t x \t \t \t J(x)\n"); */
  for (x=0.0; x<10.1; x+=0.1)
    {
      alpha = 1;
      n = 10;
      oldx = x;
      besj_(&x, &alpha, &n, &y, &nz);
      x = oldx;
      /*      printf("\t%f \t \t %f\n",x,y ); */
    }

  for (x=0.0; x<10.1; x+=0.1)
    {
      alpha = 1;
      n = 10;
      oldx = x;
      besjp(&x, &alpha, &n, &y, &nz);
      x = oldx;
      /* printf("%f \t \t %f\n",x,y ); */
    }

  for (x=0.1; x<10.1; x+=0.1)
    {
      n = 10;
      oldx = x;
      fnu = 1.0;
      besy_(&x, &fnu, &n, &y);
      x = oldx;
      /* printf("%f \t \t %f\n",x,y ); */
    }

  for (x=0.1; x<10.1; x+=0.1)
    {
      n = 10;
      oldx = x;
      fnu = 1.0;
      besyp(&x, &fnu, &n, &y);
      x = oldx;
      /*      printf("%f \t \t %f\n",x,y ); */
    }

  for (x=0.1; x<10.1; x+=0.1)
    {
      n = 10;
      oldx = x;
      alpha = 1.0;
      kode = 2;
      besi_(&x, &alpha, &kode, &n, &y, &nz);
      x = oldx;
      /*      printf("%f \t \t %f\n",x,y ); */
    }

  for (x=0.1; x<10.1; x+=0.1)
    {
      n = 10;
      oldx = x;
      alpha = 1.0;
      kode = 2;
      besip(&x, &alpha, &kode, &n, &y, &nz);
      x = oldx;
      /* printf("%f \t \t %f\n",x,y[0] ); */
    }

  for (x=0.1; x<10.1; x+=0.1)
    {
      n = 10;
      oldx = x;
      fnu = 1.0;
      kode = 2;
      besk_(&x, &fnu, &kode, &n, &y, &nz);
      x = oldx;
      /* printf("%f \t \t %f\n",x,y[0] ); */
    }

  for (x=0.1; x<10.1; x+=0.1)
    {
      n = 10;
      oldx = x;
      fnu = 1.0;
      kode = 2;
      beskp(&x, &fnu, &kode, &n, &y, &nz);
      x = oldx;
      /* printf("%f \t \t %f\n",x,y[0] ); */
    }

  fnu = 1.0;
  beta = 5.0;
  k0 = 2*M_PI/1.28;
  index = 1.3;
  wrb0 = 2.5;
  wrb1 = 3.0;
  transfer(beta, fnu, index, wrb0, wrb1, k0, ftm);
  for (i=0; i<4; i++)
    {
      for(j=0; j<4; j++)
	{
	  /*  printf("ftm[%d][%d]\t %f\n",i,j,ftm[i][j] ); */
	}
    }


  for (z__.r=2.8; z__.r<12.0; z__.r+=0.1)
    {
      oldzr = z__.r;
      n = 20;
      kode = 2;
      fnu = 1.0;
      /*      cbesj_(&z__, &fnu, &kode, &n, &cy, &nz, &ierr);
      z__.r = oldzr;
      z__.i = 1.0;
      printf("%f \t %f \t %f \t %f \t %d \n",z__.r, z__.i, cy.r, cy.i, ierr); */
    }

  nbound = 4;
  beta = 5.0;
  k0 = 2*M_PI/1.28;
  Aindex[0] = 1.454;
  Aindex[1] = 1.454;
  Aindex[2] = 1.446;
  Aindex[3] = 1.448;
  Aindex[4] = 1.448;
  Awrb[0] = 0.0;
  Awrb[1] = 2.5;  
  Awrb[2] = 3.0;
  Awrb[3] = 3.3;
  Awrb[4] = 3.5;
  /* product(nbound, beta, fnu, k0, Aindex, Awrb, M_product); */
  for (i=0;i<4;i++)
    {
      for (j=0;j<4;j++)
	{
	  /* printf("M_product[%d][%d] = %f \n",i,j,M_product[i][j]); */
	}
    }
  
  printf("k0 %f \n", k0);
  nbound = 4;
  fnu = 1.0;
  wbeta[0] = Aindex[2]*k0*1.001;
  wbeta[0] = 7.12;
  printf("wbeta = %f \n",wbeta[0]);
  printf("wbeta = %f \n",Aindex[1]*k0*0.999);
  for (i=1; i<100; i++)
    {
      wbeta[i] = wbeta[i-1] + 0.00001;
    }
      
  bigdet(nbound, wbeta, fnu, k0, Aindex, Awrb, &bigdetc);
  for (i=0; i<100; i++)
    {
      printf("wbeta[%d] = %f bigdetc[%d] = %lf \n",i,wbeta[i],i,bigdetc[i]/1e-9);
    }

}
