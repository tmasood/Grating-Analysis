/* ctrti2.f -- translated by f2c (version 19960717).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static complex c_b1 = {(float)1.,(float)0.};
static integer c__1 = 1;

/* Subroutine */ int ctrti2_(char *uplo, char *diag, integer *n, complex *a,
			     integer *lda, integer *info, ftnlen uplo_len,
			     ftnlen diag_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    complex q__1;

    /* Builtin functions */
    void c_div();

    /* Local variables */
    static integer j;
    extern /* Subroutine */ int cscal_();
    extern logical lsame_();
    static logical upper;
    extern /* Subroutine */ int ctrmv_(), xerbla_();
    static logical nounit;
    static complex ajj;


/*  -- LAPACK routine (version 3.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  CTRTI2 computes the inverse of a complex upper or lower triangular */
/*  matrix. */

/*  This is the Level 2 BLAS version of the algorithm. */

/*  Arguments */
/*  ========= */

/*  UPLO    (input) CHARACTER*1 */
/*          Specifies whether the matrix A is upper or lower triangular. 
*/
/*          = 'U':  Upper triangular */
/*          = 'L':  Lower triangular */

/*  DIAG    (input) CHARACTER*1 */
/*          Specifies whether or not the matrix A is unit triangular. */
/*          = 'N':  Non-unit triangular */
/*          = 'U':  Unit triangular */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  A       (input/output) COMPLEX array, dimension (LDA,N) */
/*          On entry, the triangular matrix A.  If UPLO = 'U', the */
/*          leading n by n upper triangular part of the array A contains 
*/
/*          the upper triangular matrix, and the strictly lower */
/*          triangular part of A is not referenced.  If UPLO = 'L', the */
/*          leading n by n lower triangular part of the array A contains 
*/
/*          the lower triangular matrix, and the strictly upper */
/*          triangular part of A is not referenced.  If DIAG = 'U', the */
/*          diagonal elements of A are also not referenced and are */
/*          assumed to be 1. */

/*          On exit, the (triangular) inverse of the original matrix, in 
*/
/*          the same storage format. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= max(1,N). */

/*  INFO    (output) INTEGER */
/*          = 0: successful exit */
/*          < 0: if INFO = -k, the k-th argument had an illegal value */

/*  ===================================================================== 
*/

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *info = 0;
    upper = lsame_(uplo, "U", 1L, 1L);
    nounit = lsame_(diag, "N", 1L, 1L);
    if (! upper && ! lsame_(uplo, "L", 1L, 1L)) {
	*info = -1;
    } else if (! nounit && ! lsame_(diag, "U", 1L, 1L)) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CTRTI2", &i__1, 6L);
	return 0;
    }

    if (upper) {

/*        Compute inverse of upper triangular matrix. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (nounit) {
		i__2 = j + j * a_dim1;
		c_div(&q__1, &c_b1, &a[j + j * a_dim1]);
		a[i__2].r = q__1.r, a[i__2].i = q__1.i;
		i__2 = j + j * a_dim1;
		q__1.r = -a[i__2].r, q__1.i = -a[i__2].i;
		ajj.r = q__1.r, ajj.i = q__1.i;
	    } else {
		q__1.r = (float)-1., q__1.i = (float)0.;
		ajj.r = q__1.r, ajj.i = q__1.i;
	    }

/*           Compute elements 1:j-1 of j-th column. */

	    i__2 = j - 1;
	    ctrmv_("Upper", "No transpose", diag, &i__2, &a[a_offset], lda, &
		    a[j * a_dim1 + 1], &c__1, 5L, 12L, 1L);
	    i__2 = j - 1;
	    cscal_(&i__2, &ajj, &a[j * a_dim1 + 1], &c__1);
/* L10: */
	}
    } else {

/*        Compute inverse of lower triangular matrix. */

	for (j = *n; j >= 1; --j) {
	    if (nounit) {
		i__1 = j + j * a_dim1;
		c_div(&q__1, &c_b1, &a[j + j * a_dim1]);
		a[i__1].r = q__1.r, a[i__1].i = q__1.i;
		i__1 = j + j * a_dim1;
		q__1.r = -a[i__1].r, q__1.i = -a[i__1].i;
		ajj.r = q__1.r, ajj.i = q__1.i;
	    } else {
		q__1.r = (float)-1., q__1.i = (float)0.;
		ajj.r = q__1.r, ajj.i = q__1.i;
	    }
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

		i__1 = *n - j;
		ctrmv_("Lower", "No transpose", diag, &i__1, &a[j + 1 + (j + 
			1) * a_dim1], lda, &a[j + 1 + j * a_dim1], &c__1, 5L, 
			12L, 1L);
		i__1 = *n - j;
		cscal_(&i__1, &ajj, &a[j + 1 + j * a_dim1], &c__1);
	    }
/* L20: */
	}
    }

    return 0;

/*     End of CTRTI2 */

} /* ctrti2_ */

