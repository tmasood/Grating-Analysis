#include "f2c.h"

/* Subroutine */ int zptcon_(integer *n, doublereal *d, doublecomplex *e, 
	doublereal *anorm, doublereal *rcond, doublereal *rwork, integer *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZPTCON computes the reciprocal of the condition number (in the   
    1-norm) of a complex Hermitian positive definite tridiagonal matrix   
    using the factorization A = L*D*L**H or A = U**H*D*U computed by   
    ZPTTRF.   

    Norm(inv(A)) is computed by a direct method, and the reciprocal of   
    the condition number is computed as   
                     RCOND = 1 / (ANORM * norm(inv(A))).   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the diagonal matrix D from the   
            factorization of A, as computed by ZPTTRF.   

    E       (input) COMPLEX*16 array, dimension (N-1)   
            The (n-1) off-diagonal elements of the unit bidiagonal factor 
  
            U or L from the factorization of A, as computed by ZPTTRF.   

    ANORM   (input) DOUBLE PRECISION   
            The 1-norm of the original matrix A.   

    RCOND   (output) DOUBLE PRECISION   
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the   
            1-norm of inv(A) computed in this routine.   

    RWORK   (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The method used is described in Nicholas J. Higham, "Efficient   
    Algorithms for Computing the Condition Number of a Tridiagonal   
    Matrix", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986.   

    ===================================================================== 
  


       Test the input arguments.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    /* Builtin functions */
    double z_abs(doublecomplex *);
    /* Local variables */
    static integer i, ix;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal ainvnm;



#define RWORK(I) rwork[(I)-1]
#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*anorm < 0.) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZPTCON", &i__1);
	return 0;
    }

/*     Quick return if possible */

    *rcond = 0.;
    if (*n == 0) {
	*rcond = 1.;
	return 0;
    } else if (*anorm == 0.) {
	return 0;
    }

/*     Check that D(1:N) is positive. */

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (D(i) <= 0.) {
	    return 0;
	}
/* L10: */
    }

/*     Solve M(A) * x = e, where M(A) = (m(i,j)) is given by   

          m(i,j) =  abs(A(i,j)), i = j,   
          m(i,j) = -abs(A(i,j)), i .ne. j,   

       and e = [ 1, 1, ..., 1 ]'.  Note M(A) = M(L)*D*M(L)'.   

       Solve M(L) * x = e. */

    RWORK(1) = 1.;
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	RWORK(i) = RWORK(i - 1) * z_abs(&E(i - 1)) + 1.;
/* L20: */
    }

/*     Solve D * M(L)' * x = b. */

    RWORK(*n) /= D(*n);
    for (i = *n - 1; i >= 1; --i) {
	RWORK(i) = RWORK(i) / D(i) + RWORK(i + 1) * z_abs(&E(i));
/* L30: */
    }

/*     Compute AINVNM = max(x(i)), 1<=i<=n. */

    ix = idamax_(n, &RWORK(1), &c__1);
    ainvnm = (d__1 = RWORK(ix), abs(d__1));

/*     Compute the reciprocal condition number. */

    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }

    return 0;

/*     End of ZPTCON */

} /* zptcon_ */

