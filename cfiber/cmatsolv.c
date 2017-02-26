#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"

extern int cgetrf(integer*, integer*, complex*, integer*,
		  integer*, integer*);
extern int cgetri(integer*, complex*, integer*, integer*, complex*,
		  integer*, integer*);

#define NCA  3

int cmatsolv(complex matA[][NCA], complex *matB, complex *matC)
{
  int i,j,k;
  integer n; /* order of the matrix matA, num rows and num cols */
  integer lda; /* leading dimension of array matA */
  integer ipiv[NCA]; /* pivot matrix */
  integer lwork;
  integer info;

  complex work[NCA];
  complex p1, p2, p;
  complex czero;

  czero.r = 0.0;
  czero.i = 0.0;

  lwork = NCA;
  n = NCA;
  lda = NCA;
  k = 0;

  cgetrf(&n, &n, matA, &lda, ipiv, &info);
  cgetri(&n, matA, &lda, ipiv, work, &lwork, &info);

  for (i=0; i<NCA; i++)
    {
      matC[i] = czero;
      for (j=0; j<NCA; j++)
	{
	  p1 = c_prod(matA[i][j], matB[j]);
	  p = c_add(p1, matC[i]);	  
	  matC[i] = p;
	}
    }
  return (0);
}
