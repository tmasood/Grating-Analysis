//
//  Copyright 2015 Quantum Designs LLC, Taha Masood, Johannes Tausch
//  and Jerome Butler
//
//  Permission to use, copy, and distribute this software and its
//  documentation for any purpose with or without fee is hereby granted,
//  provided that the above copyright notice appear in all copies and
//  that both the copyright notice and this permission notice appear
//  in supporting documentation.
//
//  This software is provided "as is" without express or implied warranty
//  to the extent permitted by applicable law.
//
#include <iostream>
#include <complex>
#include <cstdlib>

using namespace std;

extern "C" {
// SVD. Note that the routine returns v**h, not v
void zgesvd_(char *jobu, char *jobvt, int* m, int* n,
	     complex<double> *a, int* lda, double *s,
	     complex<double> *u, int* ldu, complex<double> *vt,
	     int* ldvt, complex<double> *work, int* lwork,
	     double *rwork, int* info);
}

// allocate and calculate the nullspace of the 4*order+2 square
// matrix T via SVD

// tol is the tolerance

complex<double> **getnullspace(int order, int *dimnp, complex<double> *t,
			       double tol)
{
  int i, j, k, dimn;

  complex<double> **mat;

  // don't calc left sing vects, write right ones into T
  char jobu='N', jobvt='O';  

  // number of rows of a
  int n = (4*order+2);
  // number of columns of a
  int m = n;
  // matrices u, vt                
  complex<double> *u=NULL, *vt=NULL;
  // leading dimensions of a, u, vt
  int lda, ldu, ldvt;       

  // singular values of a, sorted in descending order
  double *sing;

  // workspace
  complex<double> *work;
  // workspace   
  double *rwork; 
  // size of work, rwork           
  int lwork=10*n;
  // return flag, if zero successful exit
  int info;

  lda = ldu = ldvt = n;
  rwork = new double[lwork];
  work = new complex<double>[lwork];
  sing = new double[n];

  zgesvd_( &jobu, &jobvt, &m, &n, t, &lda, sing, u, &ldu, vt,
	   &ldvt, work, &lwork, rwork, &info);

  if (info)
    {
      cout << "getnullspace(): SVD bombed!" << endl;
      exit(1);
    }

  for (dimn=0, i=n-1; i>=0; i--, dimn++)
    {
      if (sing[i] > tol) break;
    }

  if (!dimn)
    {
      cout << "\n getnullspace(): trivial nullspace \n" << endl;
      exit(1);
    }

  mat = new complex<double>*[dimn];
  for (i = 0; i < dimn; i++)
    {
      k = n-i-1;
      mat[i] = new complex<double>[n];
      for (j = 0; j < n; j++, k+=n)
	{
	  mat[i][j] = conj(t[k]);
	}
  }

  *dimnp = dimn;
  
  return mat;

} // getnullspace
