#include "f2c.h"

/* Subroutine */ int dpbsv_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb, 
	integer *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPBSV computes the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N symmetric positive definite band matrix and X   
    and B are N-by-NRHS matrices.   

    The Cholesky decomposition is used to factor A as   
       A = U**T * U,  if UPLO = 'U', or   
       A = L * L**T,  if UPLO = 'L',   
    where U is an upper triangular band matrix, and L is a lower   
    triangular band matrix, with the same number of superdiagonals or   
    subdiagonals as A.  The factored form of A is then used to solve the 
  
    system of equations A * X = B.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)   
            On entry, the upper or lower triangle of the symmetric band   
            matrix A, stored in the first KD+1 rows of the array.  The   
            j-th column of A is stored in the j-th column of the array AB 
  
            as follows:   
            if UPLO = 'U', AB(KD+1+i-j,j) = A(i,j) for max(1,j-KD)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(N,j+KD). 
  
            See below for further details.   

            On exit, if INFO = 0, the triangular factor U or L from the   
            Cholesky factorization A = U**T*U or A = L*L**T of the band   
            matrix A, in the same storage format as A.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i of A is not   
                  positive definite, so the factorization could not be   
                  completed, and the solution has not been computed.   

    Further Details   
    ===============   

    The band storage scheme is illustrated by the following example, when 
  
    N = 6, KD = 2, and UPLO = 'U':   

    On entry:                       On exit:   

        *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46   
        *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56   
       a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66   

    Similarly, if UPLO = 'L' the format of A is as follows:   

    On entry:                       On exit:   

       a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66   
       a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *   
       a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *   

    Array elements marked * are not used by the routine.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer ab_dim1, ab_offset, b_dim1, b_offset, i__1;
    /* Local variables */
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int xerbla_(char *, integer *), dpbtrf_(
	    char *, integer *, integer *, doublereal *, integer *, integer *), dpbtrs_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *);



#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kd < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*ldab < *kd + 1) {
	*info = -6;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPBSV ", &i__1);
	return 0;
    }

/*     Compute the Cholesky factorization A = U'*U or A = L*L'. */

    dpbtrf_(uplo, n, kd, &AB(1,1), ldab, info);
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

	dpbtrs_(uplo, n, kd, nrhs, &AB(1,1), ldab, &B(1,1), ldb, 
		info);

    }
    return 0;

/*     End of DPBSV */

} /* dpbsv_ */
