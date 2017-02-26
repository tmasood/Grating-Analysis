#include <stdio.h>
#include <math.h>
#include <string.h>
#include "f2c.h"
#include "cmplx.h"

#define PI M_PI
#define SIZE 4
#define MAXSIZE 100
#define DIMENSION 200

extern int cbigdet(int, complex*, real, complex, complex*, real*, complex*);
extern int cfirstcoeff(int, complex, real, complex, complex*, real*, complex*);
extern int cmatmul(complex *, complex *, complex *);
extern int cmatmul2(complex *, complex *, complex *);
extern int ctransfer(complex, real, complex, complex, real, complex, complex *);
extern int assign_region(int, real *, real *, int *);
extern int cfieldez(complex, complex, real, real *, int *, complex *,
            int, real *, complex *, complex *);
extern int cfieldhz(complex, complex, real, real *, int *, complex *,
            int, real *, complex *, complex *);
extern int cfielder(complex, complex, real, real *, int *, complex *,
            int, real *, complex *, complex *);
extern int cfieldhr(complex, complex, real, real *, int *, complex *,
            int, real *, complex *, complex *, complex *);
extern int cfieldephi(complex, complex, real, real *, int *, complex *,
              int, real *, complex *, complex *);
extern int cfieldhphi(complex, complex, real, real *, int *, complex *,
              int, real *, complex *, complex *, complex *);


int main(int argc, char *argv[])
{

  FILE *fp;
  FILE *fpi;
  FILE *fpimag;

  real fnu;
  complex k0;

  complex index[20];
  real wrb[20];
  complex wbeta[MAXSIZE];
  complex M_x[SIZE-1];
  complex MatP[SIZE][SIZE];
  complex M_Prod2[SIZE];
  complex step;
  complex betaR;
  complex wbetah, wbetal;
  real Zimp;
  complex v_max;
  complex lambda;
  complex p1, p2, p;
  complex czero;
  complex cone;

  int i,j;
  int nbound;
  int err;
  int kount;
  int wregion[DIMENSION];
  int qtype;
  int rootfound;
  int instatus;
  char infile[50];

  complex bigdetc[MAXSIZE];
  complex qtrans[MAXSIZE];
  complex coefArr[SIZE];
  real TypeWave[MAXSIZE];
  complex coefMat[SIZE][MAXSIZE];
  real xF[DIMENSION];
  complex Ez[DIMENSION], Hz[DIMENSION];
  complex Ephi[DIMENSION], Hphi[DIMENSION];
  complex Er[DIMENSION], Hr[DIMENSION];
  complex poyn[DIMENSION];
  int state, state1, state2, state3, state4;
  integer two;

  two = 2;
  czero.r = 0.0;
  czero.i = 0.0;

  cone.r = 1.0;
  cone.i = 0.0;

  *infile = '\0';

  if (argc > 1)
    {
      strcat(infile,"dat/");
      strcat(infile,argv[1]);
      strcat(infile,".data");
    }
  else
    {
      printf("Input file not specified -- Try again \n");
      exit(1);
    }
  i = 0;
  fp = fopen("dat/field.data","w");
  fpimag = fopen("dat/fieldimag.data","w");
  fpi = fopen(infile,"r");

  if((fp == NULL) || (fpi == NULL) || (fpimag == NULL) )
    {
      printf("Cannot open file \n");
      exit(1);
    }

  fscanf(fpi,"%f\n",&lambda.r);
  printf("lambda = %f \n",lambda.r);

  k0.r = (2*M_PI)/lambda.r;
  k0.i = 0.0;
  i = 0;
  instatus = fscanf(fpi,"%f\t%f\t%f\n",&wrb[i],&index[i].r,&index[i].i);
  while (instatus != EOF)
    {
      printf("wrb %f\t index %f  %f\n",wrb[i],index[i].r,index[i].i);
      i++;
      instatus = fscanf(fpi,"%f\t%f\t%f\n",&wrb[i],&index[i].r,&index[i].i);
    }

  nbound = i - 1;
  fnu = 1.0;
  betaR = czero;
  wbetah = czero;
  wbetal = czero;

  p1 = c_prod(index[nbound], k0);
  wbeta[0] = cscalar_prod(p1, 1.000001);

  c_div(&p, &wbeta[0], &k0);
  printf("neffi.r = %f  neffi.i = %f \n",p.r, p.i);
  v_max = index[0];
  for (i=1; i<=nbound; i++)
    {
      state = c_comparel(v_max, index[i]);
      /* v_max < index[i] */
      if (state)
	{
	  v_max = index[i];
	}
    }
  printf("nefff.r = %f  nefff.i = %f \n",v_max.r, v_max.i);

  p1 = c_prod(v_max, k0);
  p2 = cscalar_prod(p1, 0.9999999);
  p1 = c_prod(index[nbound], k0);
  p1 = cscalar_prod(p1, 1.000001);
  p = c_sub(p2, p1);
  p = cscalar_div(p, 100.0);
  step = p;
  printf("step = %f  %f \n",step.r, step.i);
  for (i=1; i<100; i++)
    {
      p = c_add(wbeta[i-1], step);
      wbeta[i] = p;
    }

  /* The determinental function near a true root must */
  /* go through zero smoothly and monotonically. */

  /* monotonically: being a mathematical function that either */
  /* never decreases or never increases as the independent */
  /* variable increases   */
  err = cbigdet(nbound, wbeta, fnu, k0, index, wrb, bigdetc);

  rootfound = 0;
  for (i=2; i<100; i++)
    {
      c_div(&p, &wbeta[i], &k0);
      printf("wbeta[%d] = %f  %f bigdetc[%d] = %f  %f \n",
	     i,wbeta[i].r, wbeta[i].i, i, bigdetc[i].r, bigdetc[i].i);

      state1 = c_comparegi(bigdetc[i], czero);
      state2 = c_compareli(bigdetc[i-1], czero);
      state3 = c_compareli(bigdetc[i], czero);
      state4 = c_comparegi(bigdetc[i-1], czero);

      if (((state1) && (state2)) || 
	  ((state3) && (state4)) )
	{
	  betaR.i = (wbeta[i].i + wbeta[i-1].i)/2;
          betaR.r = (wbeta[i].r + wbeta[i-1].r)/2;
          wbetah.r = wbeta[i].r;
          wbetal.r = wbeta[i-1].r;
          wbetah.i = wbeta[i].i;
          wbetal.i = wbeta[i-1].i;
	  
	  rootfound = 1;
	  break;
	}
      else
	{
	  rootfound = 0;
	}
    }

  printf("rootfound = %d \t %f  %f\n",rootfound, betaR.r, betaR.i);
  if (rootfound == 0)
    { 
      for (i=1; i<100; i++)
	{
	  /*	  printf("wbeta[%d] = %f bigdetc[%d] = %f neff = %f\n",
		 i,wbeta[i],i,bigdetc[i],wbeta[i]/k0); */
	  if (abs(bigdetc[i].r) > abs(bigdetc[i-1].r))
	    {
	      p1 = c_add(wbeta[i], wbeta[i-1]);
	      p = cscalar_div(p1, 2.0);
	      betaR = p;
	      wbetah = wbeta[i];
	      wbetal = wbeta[i-1];
	      rootfound = 1;
	      break;
	    }
	  else
	    {
	      rootfound = 0;
	    }
	}
    }

  if (rootfound == 0)
    {
      printf("Root not found \n");
      exit(1);
    }

  for (j=0; j<2; j++)
    {
      rootfound = 0;
      p1 = c_sub(wbetah, wbetal);
      p = cscalar_div(p1, 100.0);
      step = p;

      wbeta[0] = c_add(wbetal, step);
      for (i=1; i<100; i++)
	{
	  /* wbeta[i] = c_add(wbeta[i-1], step); */
	  wbeta[i].r = wbeta[i-1].r + 7*step.r;
	  wbeta[i].i = wbeta[i-1].i + 90*step.i;
	}
      err = cbigdet(nbound, wbeta, fnu, k0, index, wrb, bigdetc);
      for (i=1; i<100; i++)
	{
	  c_div(&p, &wbeta[i], &k0);
	  printf("wbeta[%d] = %f  %f bigdetc[%d] = %f  %f \n",
	     i,wbeta[i].r, wbeta[i].i, i, bigdetc[i].r, bigdetc[i].i);

	  state1 = c_compareg(bigdetc[i], czero);
	  state2 = c_comparel(bigdetc[i-1], czero);
	  state3 = c_comparel(bigdetc[i], czero);
	  state4 = c_compareg(bigdetc[i-1], czero);
	  if (((state1) && (state2)) ||
	      ((state3) && (state4)) )
	    {
	      p1 = c_add(wbeta[i], wbeta[i-1]);
	      p = cscalar_div(p1, 2.0);
	      betaR = p;
	      wbetah = wbeta[i];
	      wbetal = wbeta[i-1];
	      rootfound = 1;
	      break;
	    }
	}
      
      if (rootfound == 0)
	{ 
	  for (i=1; i<100; i++)
	    {
	      /* printf("wbeta[%d] = %.25f bigdetc[%d] = %.25f neff = %f\n",
		     i,wbeta[i],i,bigdetc[i],wbeta[i]/k0); */
	      if (abs(bigdetc[i].r) > abs(bigdetc[i-1].r))
		{
		  p1 = c_add(wbeta[i], wbeta[i-1]);
		  p = cscalar_div(p1, 2.0);
		  betaR = p;
		  wbetah = wbeta[i];
		  wbetal = wbeta[i-1];
		  rootfound = 1;
		  break;
		}
	      else
		{
		  rootfound = 0;
		}
	    }
	}
    }

  c_div(&p, &betaR, &k0);
  printf("BetaR = %f %f neff = %f   %f \n",betaR.r, betaR.i, betaR.r/k0.r, k0.i);
  cfirstcoeff(nbound, betaR, fnu, k0, index, wrb, M_x);

  p1 = c_prod(index[0], k0);
  pow_ci(&p2, &p1, &two);
  pow_ci(&p1, &betaR, &two);
  p = c_sub(p2, p1);
  state = c_comparege(p, czero);
  if (state)
    {
      coefArr[0] = cone;
      coefArr[1] = czero;
      coefArr[2] = M_x[0];
      coefArr[3] = czero;
      TypeWave[0] = 1;
    }
  else
    {
      coefArr[0] = czero;
      coefArr[1] = cone;
      coefArr[2] = czero;
      coefArr[3] = M_x[0];
      TypeWave[0] = -1;
    }
  for (i=0; i<4; i++)
    {
      for (j=0; j<MAXSIZE; j++)
	{
	  coefMat[i][j] = czero;
	}
    }
  for (i=0; i<4; i++)
    {
      coefMat[i][0] = coefArr[i];
    }

  coefMat[0][nbound] = M_x[1];
  coefMat[1][nbound] = czero;
  coefMat[2][nbound] = M_x[2];
  coefMat[3][nbound] = czero;

  kount = 0;
  while (kount < (nbound - 1))
    {
      qtype = ctransfer(betaR, fnu, index[kount+1], index[kount],
			wrb[kount+1], k0, MatP);
      for (i=0; i<4; i++)
	{
	  for (j=0; j<4; j++)
	    {
	      printf("MatP = %f  %f \n",MatP[i][j].r, MatP[i][j].i);
	    }
	}
      cmatmul2(MatP, coefArr, M_Prod2);

      for (i=0; i<4; i++)
	{
	  coefArr[i] = M_Prod2[i];
	  printf("M_Prod2 = %f   %f \n",M_Prod2[i].r,M_Prod2[i].i);
	  coefMat[i][kount+1] = coefArr[i];
	}

      TypeWave[kount + 1] = qtype;
      kount++;
    }

  TypeWave[nbound] = -1;
  
  i = 0;
  while (i < (nbound + 1))
    {
      if (TypeWave[i] > 0)
	{
	  p1 = c_prod(k0, index[i]);
	  pow_ci(&p2, &p1, &two);
	  pow_ci(&p1, &betaR, &two);
	  p1 = c_sub(p2, p1);
	  c_sqrt(&p, &p1);
	  qtrans[i] = p;
	}
      else
	{
	  p1 = c_prod(k0, index[i]);
	  pow_ci(&p2, &p1, &two);
	  pow_ci(&p1, &betaR, &two);
	  p1 = c_sub(p1, p2);
	  c_sqrt(&p, &p1);
	  qtrans[i] = p;
	}
      i++;
    }
  
  xF[0] = 0.06;
  for (i=1; i<200; i++)
    {
      xF[i] = xF[i - 1] + 0.06;
    }

  assign_region(nbound, xF, wrb, wregion);
  cfieldez(betaR, k0, fnu, xF, wregion, qtrans,
	  nbound, TypeWave, coefMat, Ez);
  cfieldhz(betaR, k0, fnu, xF, wregion, qtrans,
	  nbound, TypeWave, coefMat, Hz);
  cfieldephi(betaR, k0, fnu, xF, wregion, qtrans,
          nbound, TypeWave, coefMat, Ephi);
  cfieldhphi(betaR, k0, fnu, xF, wregion, qtrans,
          nbound, TypeWave, coefMat, index, Hphi);
  cfielder(betaR, k0, fnu, xF, wregion, qtrans,
          nbound, TypeWave, coefMat, Er);
  cfieldhr(betaR, k0, fnu, xF, wregion, qtrans,
          nbound, TypeWave, coefMat, index, Hr);

  Zimp = 120*PI;
  for (i=0; i<200; i++)
    {
      p1 = c_prod(Er[i], Hphi[i]);
      p2 = c_prod(Ephi[i], Hr[i]);
      p = c_sub(p1, p2);
      p = cscalar_prod(p, 0.5);
      p = cscalar_div(p, Zimp);
      poyn[i] = p;
    }

  for (i=1; i<200; i++)
    {
      fprintf(fp,"%f \t %f \t %f \t %f \t %f \t %f \t %f \n",
	      xF[i],Ez[i].r, Hz[i].r,
	      Ephi[i].r, Hphi[i].r, Er[i].r, Hr[i].r);

      fprintf(fpimag,"%f \t %f \t %f \t %f \t %f \t %f \t %f \n",
	      xF[i],Ez[i].i, Hz[i].i,
	      Ephi[i].i, Hphi[i].i, Er[i].i, Hr[i].i);
    }

  fclose(fp);
  fclose(fpimag);
  fclose(fpi);

  return(0);
}
