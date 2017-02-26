#include<stdio.h>
#include<math.h>
#include "f2c.h"

extern int cbesj(complex*, real*, integer*, integer*, complex*, integer*, integer*);
extern int cbesyp(complex*, real*, integer*, integer, real*);
extern int cbesjp(real*, real*, int*, real*, int*);
extern int cbesy(real*, real*, int*, real*);
extern int cbesip(real*, real*, int*, int*, real*, int*);
extern int cbesk(real*, real*, int*, int*, real*, int*);
extern int cbesi(real*, real*, int*, int*, real*, int*);
extern int cbeskp(real*, real*, int*, int*, real*, int*);

#define PI M_PI
#define SIZE 4

int transfer(real beta, real nu, real nref, real r1, real r2,
	     real k0, real ftm[][SIZE])
{
  real lay[SIZE][SIZE];
  real kappa, gamma;
  real alpha1;
  real alpha2;
  real rho;
  real sigma;
  real z1;
  real z2;
  real del1, del2, del3, del4;
  real alpha;
  real fnu;
  real del1bjy[10], del1bjpy[10];
  real del1byy[10], del1bypy[10];
  real del2bjy[10], del2bjpy[10];
  real del2byy[10], del2bypy[10];
  real del3bjy1[10], del3bjy2[10];
  real del3byy1[10], del3byy2[10];
  real del4bjpy1[10], del4bjpy2[10];
  real del4bypy1[10], del4bypy2[10];

  real del1bipy[10], del1bky[10], del1biy[10];
  real del1bkpy[10], del2bipy[10], del2bky[10],del2bkpy[10];
  real del2biy[10], del3biy2[10], del3bky1[10];
  real del3biy1[10], del3bky2[10], del4bipy2[10];
  real del4bkpy1[10], del4bipy1[10], del4bkpy2[10];

  int nz, n;
  int kode;
  int i, j;

  kode = 1;
  n = 10;
  alpha = nu;
  fnu = nu;

  if (nref*k0 >= beta)
    {
      kappa = sqrt(pow((nref*k0),2) - pow(beta,2));
      alpha1 = (beta * nu) / (pow(kappa,2) * r1);
      alpha2 = (beta * nu) / (pow(kappa,2) * r2);
      rho = k0/kappa;
      sigma = (k0 * pow(nref,2))/kappa;
      z1 = kappa * r1;
      z2 = kappa * r2;
      besj(&z2, &alpha, &n, del1bjy, &nz);
      besyp(&z1, &fnu, &n, del1bypy);
      besjp(&z1, &alpha, &n, del1bjpy, &nz);
      besy(&z2, &fnu, &n, del1byy);
      del1 = (del1bjy[0] * del1bypy[0]) - (del1bjpy[0] * del1byy[0]);
      
      besj(&z1, &alpha, &n, del2bjy, &nz);
      besyp(&z2, &fnu, &n, del2bypy);
      besjp(&z2, &alpha, &n, del2bjpy, &nz);
      besy(&z1, &fnu, &n, del2byy);
      del2 = (del2bjy[0] * del2bypy[0]) - (del2bjpy[0] * del2byy[0]);

      besj(&z2, &alpha, &n, del3bjy2, &nz);
      besy(&z1, &fnu, &n, del3byy1);
      besj(&z1, &alpha, &n, del3bjy1, &nz);
      besy(&z2, &fnu, &n, del3byy2);
      del3 = (del3bjy2[0] * del3byy1[0]) - (del3bjy1[0] * del3byy2[0]);

      besjp(&z2, &alpha, &n, del4bjpy2, &nz);
      besyp(&z1, &fnu, &n, del4bypy1);
      besjp(&z1, &alpha, &n, del4bjpy1, &nz);
      besyp(&z2, &fnu, &n, del4bypy2);
      del4 = (del4bjpy2[0] * del4bypy1[0]) - (del4bjpy1[0] * del4bypy2[0]);


      lay[0][0] = del1;
      lay[0][1] = (alpha1*del3)/sigma;
      lay[0][2] = 0.0;
      lay[0][3] = -del3/sigma;
      
      lay[1][0] = (alpha1*del3)/rho;
      lay[1][1] = del1;
      lay[1][2] = -del3/rho;
      lay[1][3] = 0.0;

      lay[2][0] = (alpha2*del1) - (alpha1*del2);
      lay[2][1] = ((alpha1*alpha2*del3)/sigma) + (rho*del4);
      lay[2][2] = del2;
      lay[2][3] = -(alpha2*del3)/sigma;
      
      lay[3][0] = ((alpha1*alpha2*del3)/rho) + (sigma*del4);
      lay[3][1] = (alpha2*del1) - (alpha1*del2);
      lay[3][2] = -(alpha2*del3)/rho;
      lay[3][3] = del2;

      for (i=0;i<4;i++)
	{
	  for (j=0;j<4;j++)
	    {
	      ftm[i][j] = lay[i][j] * ((PI * z1)/2.0);
	    }
	}
      return (1);
    }
  
  else
    {
      gamma = sqrt(pow(beta,2) - pow((nref*k0),2));
      alpha1 = (beta * nu) / (pow(gamma,2) * r1);
      alpha2 = (beta * nu) / (pow(gamma,2) * r2);
      rho = k0/gamma;
      sigma = (k0 * pow(nref,2))/gamma;
      z1 = gamma * r1;
      z2 = gamma * r2;
      besip(&z1, &alpha, &kode, &n, del1bipy, &nz);
      besk(&z2, &fnu, &kode, &n, del1bky, &nz);
      besi(&z2, &alpha, &kode, &n, del1biy, &nz);
      beskp(&z1, &fnu, &kode, &n, del1bkpy, &nz);
      del1 = (del1bipy[0] * del1bky[0]) - (del1biy[0] * del1bkpy[0]);
      
      besip(&z2, &alpha, &kode, &n, del2bipy, &nz);
      besk(&z1, &fnu, &kode, &n, del2bky, &nz);
      besi(&z1, &alpha, &kode, &n, del2biy, &nz);
      beskp(&z2, &fnu, &kode, &n, del2bkpy, &nz);
      del2 = (del2bipy[0] * del2bky[0]) - (del2biy[0] * del2bkpy[0]);

      besi(&z1, &alpha, &kode, &n, del3biy2, &nz);
      besk(&z2, &fnu, &kode, &n, del3bky1, &nz);
      besi(&z2, &alpha, &kode, &n, del3biy1, &nz);
      besk(&z1, &fnu, &kode, &n, del3bky2, &nz);
      del3 = (del3biy2[0] * del3bky1[0]) - (del3biy1[0] * del3bky2[0]);

      besip(&z1, &alpha, &kode, &n, del4bipy2, &nz);
      beskp(&z2, &fnu, &kode, &n, del4bkpy1, &nz);
      besip(&z2, &alpha, &kode, &n, del4bipy1, &nz);
      beskp(&z1, &fnu, &kode, &n, del4bkpy2, &nz);
      del4 = (del4bipy2[0] * del4bkpy1[0]) - (del4bipy1[0] * del4bkpy2[0]);

      lay[0][0] = del1;
      lay[0][1] = (alpha1*del3)/sigma;
      lay[0][2] = 0.0;
      lay[0][3] = del3/sigma;
      
      lay[1][0] = (alpha1*del3)/rho;
      lay[1][1] = del1;
      lay[1][2] = del3/rho;
      lay[1][3] = 0.0;

      lay[2][0] = -(alpha2*del1) + (alpha1*del2);
      lay[2][1] = -((alpha1*alpha2*del3)/sigma) - (rho*del4);
      lay[2][2] = del2;
      lay[2][3] = -(alpha2*del3)/sigma;
      
      lay[3][0] = -((alpha1*alpha2*del3)/rho) - (sigma*del4);
      lay[3][1] = -(alpha2*del1) + (alpha1*del2);
      lay[3][2] = -(alpha2*del3)/rho;
      lay[3][3] = del2;

      for (i=0;i<4;i++)
	{
	  for (j=0;j<4;j++)
	    {
	      ftm[i][j] = lay[i][j] * z1;
	    }
	}
      return (-1);
    }
}
