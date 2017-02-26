#include<stdio.h>
#include<math.h>
#include "f2c.h"
 
#define PI M_PI
#define SIZE 4
#define DIMENSION 200
#define MAXSIZE 100

/* find Coefficient transfer matrix for a given layer */
int fex(real x, real y, real betaroot, real k0, real nu,
	real *xF, int *wregion, real *Qtrans, int nbound, 
	real *typewave, real coefMat[][MAXSIZE], real *Ez)
{
  real ztype[DIMENSION];
  real alpha;
  real z2;
  real fnu;
  real besjqt[10], besyqt[10];
  real beskqt[10], besiqt[10];
  real Qt;
  
  int numx;
  int numr;
  int region;
  int i;
  int nz, n;
  int kode;

  n = 10;
  kode = 1;
  alpha = 1.0;
  fnu = 1.0;
  numx = DIMENSION;
  numr = nbound + 1;
  i = 0;

  r = sqrt(pow(x,2) + pow(y,2));
  theta = atan(y/x);

  return (0);
}
