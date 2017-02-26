#include "f2c.h"

/* Subroutine */ int strti2_(char *uplo, char *diag, integer *n, real *a, 
	integer *lda, integer *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    STRTI2 computes the inverse of a real upper or lower triangular   
    matrix.   

    This is the Level 2 BLAS version of the algorithm.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the matrix A is upper or lower triangular. 
  
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    DIAG    (input) CHARACTER*1   
            Specifies whether or not the matrix A is unit triangular.   
            = 'N':  Non-unit triangular   
            = 'U':  Unit triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) REAL array, dimension (LDA,N)   
            On entry, the triangular matrix A.  If UPLO = 'U', the   
            leading n by n upper triangular part of the array A contains 
  
            the upper triangular matrix, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n by n lower triangular part of the array A contains 
  
            the lower triangular matrix, and the strictly upper   
            triangular part of A is not referenced.  If DIAG = 'U', the   
            diagonal elements of A are also not referenced and are   
            assumed to be 1.   

            On exit, the (triangular) inverse of the original matrix, in 
  
            the same storage format.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer j;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static logical upper;
    extern /* Subroutine */ int strmv_(char *, char *, char *, integer *, 
	    real *, integer *, real *, integer *), 
	    xerbla_(char *, integer *);
    static logical nounit;
    static real ajj;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    nounit = lsame_(diag, "N");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (! nounit && ! lsame_(diag, "U")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("STRTI2", &i__1);
	return 0;
    }

    if (upper) {

/*        Compute inverse of upper triangular matrix. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    if (nounit) {
		A(j,j) = 1.f / A(j,j);
		ajj = -(doublereal)A(j,j);
	    } else {
		ajj = -1.f;
	    }

/*           Compute elements 1:j-1 of j-th column. */

	    i__2 = j - 1;
	    strmv_("Upper", "No transpose", diag, &i__2, &A(1,1), lda, &
		    A(1,j), &c__1);
	    i__2 = j - 1;
	    sscal_(&i__2, &ajj, &A(1,j), &c__1);
/* L10: */
	}
    } else {

/*        Compute inverse of lower triangular matrix. */

	for (j = *n; j >= 1; --j) {
	    if (nounit) {
		A(j,j) = 1.f / A(j,j);
		ajj = -(doublereal)A(j,j);
	    } else {
		ajj = -1.f;
	    }
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

		i__1 = *n - j;
		strmv_("Lower", "No transpose", diag, &i__1, &A(j+1,j+1), lda, &A(j+1,j), &c__1);
		i__1 = *n - j;
		sscal_(&i__1, &ajj, &A(j+1,j), &c__1);
	    }
/* L20: */
	}
    }

    return 0;

/*     End of STRTI2 */

} /* strti2_ */

