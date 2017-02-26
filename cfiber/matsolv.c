#include<stdio.h>
#include<math.h>
#include "f2c.h"

extern int sgetrf(int*, int*, real*, int*, int*, int*);
extern int sgetri(int*, real*, int*, int*, real*, int*, int*);

#define NCA  3

int matsolv(real matA[][NCA], real *matB, real *matC)
{
  int i,j,k;
  int n; /* order of the matrix matA, num rows and num cols */
  int lda; /* leading dimension of array matA */
  int ipiv[NCA]; /* pivot matrix */
  int lwork;
  int info;

  real work[NCA];

  lwork = NCA;
  n = NCA;
  lda = NCA;
  k = 0;

  sgetrf(&n, &n, matA, &lda, &ipiv, &info);
  sgetri(&n, matA, &lda, &ipiv, &work, &lwork, &info);

  for (i=0; i<NCA; i++)
    {
      matC[i] = 0.0;
      for (j=0; j<NCA; j++)
	{
	  matC[i] = matC[i] + (matA[i][j] * matB[j]);
	}
    }
  return (0);
}
