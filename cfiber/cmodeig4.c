#include <stdio.h>
#include <math.h>
#include "f2c.h"

#define PI M_PI
#define SIZE 4
#define MAXSIZE 100
#define DIMENSION 200

extern int bigdet(int, real*, real, real, real*, real*, real*);
extern int firstcoeff(int, real, real, real, real*, real*, real*);
extern int matmul(real *, real *, real *);
extern int matmul2(real *, real *, real *);
extern int ctransfer(real, real, real, real, real, real, real *);
extern int assign_region(int, real *, real *, int *);
extern int fieldez(real, real, real, real *, int *, real *,
            int, real *, real *, real *);
extern int fieldhz(real, real, real, real *, int *, real *,
            int, real *, real *, real *);
extern int fieldephi(real, real, real, real *, int *, real *,
              int, real *, real *, real *);
extern int fieldhphi(real, real, real, real *, int *, real *,
              int, real *, real *, real *, real *);
extern int fielder(real, real, real, real *, int *, real *,
            int, real *, real *, real *);
extern int fieldhr(real, real, real, real *, int *, real *,
            int, real *, real *, real *, real *);

int main()
{

  FILE *fp;

  real fnu;
  real k0;
  real index[20];
  real wrb[20];
  real wbeta[MAXSIZE];
  real M_x[SIZE-1];
  real MatP[SIZE][SIZE];
  real M_Prod2[SIZE];
  real step;
  real betaR;
  real wbetah, wbetal;
  real Zimp;
  real v_max;

  int i,j;
  int nbound;
  int err;
  int kount;
  int wregion[DIMENSION];
  int qtype;

  real bigdetc[MAXSIZE];
  real qtrans[MAXSIZE];
  real coefArr[SIZE];
  real TypeWave[MAXSIZE];
  real coefMat[SIZE][MAXSIZE];
  real xF[DIMENSION];
  real Ez[DIMENSION], Hz[DIMENSION];
  real Ephi[DIMENSION], Hphi[DIMENSION];
  real Er[DIMENSION], Hr[DIMENSION];
  real poyn[DIMENSION];

  fp = fopen("field.data","w");
  if(fp == NULL)
    {
      printf("Cannot open file \n");
      exit(1);
    }

  k0 = 2*M_PI/1.5500;
  index[0] = 1.45237;
  index[1] = 1.44716;

  wrb[0] = 0.0;
  wrb[1] = 4.15;

  nbound = 1;
  fnu = 1.0;
  betaR = 0.0;
  wbetah = 0;
  wbetal = 0;
  wbeta[0] = index[nbound]*k0*1.000001;
  printf("wbetai = %f \n",wbeta[0]);
  v_max = index[0];
  for (i=1; i<=nbound; i++)
    {
      if (v_max < index[i])
	{
	  v_max = index[i];
	}
    }
  printf("wbetaf = %f \n",v_max*k0);
  step = ((v_max*k0*0.9999999) - (index[nbound]*k0*1.000001))/100;

  for (i=1; i<100; i++)
    {
      wbeta[i] = wbeta[i-1] + step;
    }
      
  err = bigdet(nbound, wbeta, fnu, k0, index, wrb, bigdetc);
  for (i=0; i<100; i++)
    {
      /*      printf("wbeta[%d] = %f bigdetc[%d] = %f \n",i,wbeta[i],i,bigdetc[i]); */
      if (((bigdetc[i] > 0) && (bigdetc[i-1] < 0)) || ((bigdetc[i] < 0) && (bigdetc[i-1] > 0)) )
	{
	  betaR = (wbeta[i] + wbeta[i-1])/2;
	  wbetah = wbeta[i];
	  wbetal = wbeta[i-1];
	}
    }

  for (j=0; j<5; j++)
    {
      step = (wbetah - wbetal)/100;
      wbeta[0] = wbetal + step;

      for (i=1; i<100; i++)
	{
	  wbeta[i] = wbeta[i-1] + step;
	}
      err = bigdet(nbound, wbeta, fnu, k0, index, wrb, bigdetc);
      for (i=0; i<100; i++)
	{
	  printf("wbeta[%d] = %.25f bigdetc[%d] = %.25f \n",i,wbeta[i],i,bigdetc[i]);
	  if (((bigdetc[i] > 0) && (bigdetc[i-1] < 0)) || ((bigdetc[i] < 0) && (bigdetc[i-1] > 0)) )
	    {
	      betaR = (wbeta[i] + wbeta[i-1])/2;
	      wbetah = wbeta[i];
	      wbetal = wbeta[i-1];
	    }
	}
    }
  printf("BetaR = %.10f neff = %.10f fnu = %f\n",betaR,betaR/k0, fnu);
  firstcoeff(nbound, betaR, fnu, k0, index, wrb, M_x);

  if (((pow(index[0]*k0,2)) - (pow(betaR,2))) >= 0)
    {
      coefArr[0] = 1;
      coefArr[1] = 0;
      coefArr[2] = M_x[0];
      coefArr[3] = 0;
      TypeWave[0] = 1;
    }
  else
    {
      coefArr[0] = 0;
      coefArr[1] = 1;
      coefArr[2] = 0;
      coefArr[3] = M_x[0];
      TypeWave[0] = -1;
    }
  for (i=0; i<4; i++)
    {
      for (j=0; j<MAXSIZE; j++)
	{
	  coefMat[i][j] = 0;
	}
    }
  for (i=0; i<4; i++)
    {
      coefMat[i][0] = coefArr[i];
    }

  coefMat[0][nbound] = M_x[1];
  coefMat[1][nbound] = 0;
  coefMat[2][nbound] = M_x[2];
  coefMat[3][nbound] = 0;

  kount = 0;
  while (kount < (nbound - 1))
    {
      qtype = ctransfer(betaR, fnu, index[kount+1],index[kount], wrb[kount+1], k0, MatP);
      matmul2(MatP, coefArr, M_Prod2);

      for (i=0; i<4; i++)
	{
	  coefArr[i] = M_Prod2[i];
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
	  qtrans[i] = sqrt(pow((k0*index[i]),2) - pow(betaR,2));
	}
      else
	{
	  qtrans[i] = sqrt(pow(betaR,2) - pow((k0*index[i]),2));
	}
      i++;
    }
  
  xF[0] = 0;
  for (i=1; i<200; i++)
    {
      xF[i] = xF[i - 1] + 0.06;
    }

  assign_region(nbound, xF, wrb, wregion);
  fieldez(betaR, k0, fnu, xF, wregion, qtrans,
	  nbound, TypeWave, coefMat, Ez);
  fieldhz(betaR, k0, fnu, xF, wregion, qtrans,
	  nbound, TypeWave, coefMat, Hz);
  fieldephi(betaR, k0, fnu, xF, wregion, qtrans,
          nbound, TypeWave, coefMat, Ephi);
  fieldhphi(betaR, k0, fnu, xF, wregion, qtrans,
          nbound, TypeWave, coefMat, index, Hphi);
  fielder(betaR, k0, fnu, xF, wregion, qtrans,
          nbound, TypeWave, coefMat, Er);
  fieldhr(betaR, k0, fnu, xF, wregion, qtrans,
          nbound, TypeWave, coefMat, index, Hr);

  for (i=0; i<200; i++)
    {
      /*      printf("%f \t %f \t %f \t %f \t %f \t %f \t %f \n",xF[i], Ez[i], Hz[i],
	     Ephi[i], Hphi[i], Er[i], Hr[i]);*/
      fprintf(fp,"%f \t %f \t %f \t %f \t %f \t %f \t %f \n",xF[i],Ez[i],Hz[i],
	      Ephi[i],Hphi[i],Er[i],Hr[i]);
    }
  Zimp = 120*PI;
  for (i=0; i<200; i++)
    {
      poyn[i] = 0.5*(Er[i]*Hphi[i] - Ephi[i]*Hr[i])/Zimp;
      /* printf("%f \t %.10f \n",xF[i], poyn[i]);  */
    }
  fclose(fp);	  
  return(0);
}
