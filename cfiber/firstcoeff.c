#include<stdio.h>
#include<math.h>
#include "f2c.h"
 
#define PI M_PI
#define SIZE 4

extern int product(int, real, real, real, real*, real*, real*);
extern int besj(real*, real*, int*, real*, int*);
extern int besjp(real*, real*, int*, real*, int*);
extern int besip(real*, real*, int*, int*, real*, int*);
extern int besk(real*, real*, int*, int*, real*, int*);
extern int besi(real*, real*, int*, int*, real*, int*);
extern int beskp(real*, real*, int*, int*, real*, int*);
extern int matmul(real*, real*, real*);
extern int matsolv(real*, real*, real*);
/* find first set of eigenmode coefficients */
int firstcoeff(int nbound, real betaR, real nu, real k0, real *index, real *wrb, real *M_x)
{
  real M_product[SIZE][SIZE];
  real temp[SIZE][SIZE];
  real EigmodeMat[SIZE - 1][SIZE - 1];
  real EigVect[SIZE - 1];
  real M0[SIZE][SIZE];
  real besjm0[10], besjpm0[10],besim0[10], besipm0[10];
  real beskcm[10], beskpcm[10];
  real gamma2;
  real gamma;
  real z2, z1;
  real kappa;
 
  int i,j;
  int n = 10;
  int nz;
  int  kode;

  kode = 1;
  
  product(nbound, betaR, nu, k0, index, wrb, M_product);
  for (i=0; i<SIZE; i++)
    {
      for (j=0; j<SIZE; j++)
	{
	  temp[i][j] = M_product[i][j];
	}
    }
  gamma2 = sqrt(pow(betaR,2) - pow((index[nbound]*k0),2));
  z2 = wrb[nbound]*gamma2;

  if ((pow((index[0]*k0),2) - pow(betaR,2)) >= 0)
    {
      kappa = sqrt(pow((index[0]*k0),2) - pow(betaR,2));
      z1 = wrb[1]*kappa;
      besj(&z1, &nu, &n, besjm0, &nz);
      besjp(&z1, &nu, &n, besjpm0, &nz);
      M0[0][0] = besjm0[0];
      M0[0][1] = 0;
      M0[0][2] = 0;
      M0[0][3] = 0;
      
      M0[1][0] = 0;
      M0[1][1] = besjm0[0];
      M0[1][2] = 0;
      M0[1][3] = 0;
      
      M0[2][0] = ((betaR*nu)/(kappa*z1))*besjm0[0];
      M0[2][1] = (k0/kappa)*besjpm0[0];
      M0[2][2] = 0;
      M0[2][3] = 0;
 
      M0[3][0] = ((k0*pow(index[0],2))/kappa)*besjpm0[0];
      M0[3][1] = (betaR*nu/(kappa*z1))*besjm0[0];
      M0[3][2] = 0;
      M0[3][3] = 0;
    }
  else
    {
      gamma = sqrt(pow(betaR,2) - pow((index[0]*k0),2));
      z1 = wrb[1]*gamma;
      besi(&z1, &nu, &kode, &n, besim0, &nz);
      besip(&z1, &nu, &kode, &n, besipm0, &nz);
      M0[0][0] = besim0[0];
      M0[0][1] = 0;
      M0[0][2] = 0;
      M0[0][3] = 0;
          
      M0[1][0] = 0;
      M0[1][1] = besim0[0];
      M0[1][2] = 0;
      M0[1][3] = 0;
      
      M0[2][0] = -((betaR*nu)/(gamma*z1))*besim0[0];
      M0[2][1] = -(k0/gamma)*besipm0[0];
      M0[2][2] = 0;
      M0[2][3] = 0;
          
      M0[3][0] = -((k0*pow(index[0],2))/gamma)*besipm0[0];
      M0[3][1] = -((betaR*nu)/(gamma*z1))*besim0[0];
      M0[3][2] = 0;
      M0[3][3] = 0;
    }
  matmul(temp, M0, M_product);
  besk(&z2, &nu, &kode, &n, beskcm, &nz);
  beskp(&z2, &nu, &kode, &n, beskpcm, &nz);

      EigmodeMat[0][0] = M_product[0][1];
      EigmodeMat[0][1] = -beskcm[0];
      EigmodeMat[0][2] = 0;
 
      EigmodeMat[1][0] = M_product[1][1];
      EigmodeMat[1][1] = 0;
      EigmodeMat[1][2] = -beskcm[0];
      
      EigmodeMat[2][0] = M_product[2][1];
      EigmodeMat[2][1] = ((betaR*nu)/(gamma2*z2))*beskcm[0];
      EigmodeMat[2][2] = (k0/gamma2)*beskpcm[0];

      /*     for (i=0;i<3;i++)
	{
	  for (j=0;j<3;j++)
	    {
	      printf("EigenModeMat[%d][%d] = %f \n",i,j,EigmodeMat[i][j]);
	    }
	}
	*/
 
      EigVect[0] = -M_product[0][0];
      EigVect[1] = -M_product[1][0];
      EigVect[2] = -M_product[2][0];

      /*	  for (j=0;j<3;j++)
	    {
	      printf("EigenVect[%d] = %f \n",j,EigVect[j]);
	    }
	    */
      matsolv(EigmodeMat, EigVect, M_x);
      return (0);
}
