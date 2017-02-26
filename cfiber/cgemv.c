/* cgemv.f -- translated by f2c (version 19960717).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* DECK CGEMV */
/* Subroutine */ int cgemv_(char *trans, integer *m, integer *n, complex *alpha,
			    complex *a, integer *lda, complex *x, integer *incx,
			    complex *beta, complex *y, integer *incy,
			    ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void r_cnjg();

    /* Local variables */
    static integer info;
    static complex temp;
    static integer lenx, leny, i__, j;
    extern logical lsame_();
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_();
    static logical noconj;

/* ***BEGIN PROLOGUE  CGEMV */
/* ***PURPOSE  Multiply a complex vector by a complex general matrix. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1B4 */
/* ***TYPE      COMPLEX (SGEMV-S, DGEMV-D, CGEMV-C) */
/* ***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA */
/* ***AUTHOR  Dongarra, J. J., (ANL) */
/*           Du Croz, J., (NAG) */
/*           Hammarling, S., (NAG) */
/*           Hanson, R. J., (SNLA) */
/* ***DESCRIPTION */

/*  CGEMV  performs one of the matrix-vector operations */

/*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or */

/*     y := alpha*conjg( A' )*x + beta*y, */

/*  where alpha and beta are scalars, x and y are vectors and A is an */
/*  m by n matrix. */

/*  Parameters */
/*  ========== */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */

/*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y. */

/*              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y. */

/*           Unchanged on exit. */

/*  M      - INTEGER. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - INTEGER. */
/*           On entry, N specifies the number of columns of the matrix A. 
*/
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - COMPLEX         . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - COMPLEX          array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. */
/*           Unchanged on exit. */

/*  LDA    - INTEGER. */
/*           On entry, LDA specifies the first dimension of A as declared 
*/
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */

/*  X      - COMPLEX          array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/*           Before entry, the incremented array X must contain the */
/*           vector x. */
/*           Unchanged on exit. */

/*  INCX   - INTEGER. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - COMPLEX         . */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - COMPLEX          array of DIMENSION at least */
/*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/*           Before entry with BETA non-zero, the incremented array Y */
/*           must contain the vector y. On exit, Y is overwritten by the 
*/
/*           updated vector y. */

/*  INCY   - INTEGER. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/* ***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and */
/*                 Hanson, R. J.  An extended set of Fortran basic linear 
*/
/*                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1, */
/*                 pp. 1-17, March 1988. */
/* ***ROUTINES CALLED  LSAME, XERBLA */
/* ***REVISION HISTORY  (YYMMDD) */
/*   861022  DATE WRITTEN */
/*   910605  Modified to meet SLATEC prologue standards.  Only comment */
/*           lines were modified.  (BKS) */
/* ***END PROLOGUE  CGEMV */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Parameters .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/* ***FIRST EXECUTABLE STATEMENT  CGEMV */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame_(trans, "N", 1L, 1L) && ! lsame_(trans, "T", 1L, 1L) && ! 
	    lsame_(trans, "C", 1L, 1L)) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("CGEMV ", &info, 6L);
	return 0;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || alpha->r == (float)0. && alpha->i == (float)0. 
	    && (beta->r == (float)1. && beta->i == (float)0.)) {
	return 0;
    }

    noconj = lsame_(trans, "T", 1L, 1L);

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
*/
/*     up the start points in  X  and  Y. */

    if (lsame_(trans, "N", 1L, 1L)) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

    if (beta->r != (float)1. || beta->i != (float)0.) {
	if (*incy == 1) {
	    if (beta->r == (float)0. && beta->i == (float)0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = i__;
		    y[i__2].r = (float)0., y[i__2].i = (float)0.;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = i__;
		    i__3 = i__;
		    q__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    q__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
		    y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (beta->r == (float)0. && beta->i == (float)0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = iy;
		    y[i__2].r = (float)0., y[i__2].i = (float)0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = iy;
		    i__3 = iy;
		    q__1.r = beta->r * y[i__3].r - beta->i * y[i__3].i, 
			    q__1.i = beta->r * y[i__3].i + beta->i * y[i__3]
			    .r;
		    y[i__2].r = q__1.r, y[i__2].i = q__1.i;
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (alpha->r == (float)0. && alpha->i == (float)0.) {
	return 0;
    }
    if (lsame_(trans, "N", 1L, 1L)) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = jx;
		if (x[i__2].r != (float)0. || x[i__2].i != (float)0.) {
		    i__2 = jx;
		    q__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    q__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
		    temp.r = q__1.r, temp.i = q__1.i;
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__;
			i__4 = i__;
			i__5 = i__ + j * a_dim1;
			q__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				q__2.i = temp.r * a[i__5].i + temp.i * a[i__5]
				.r;
			q__1.r = y[i__4].r + q__2.r, q__1.i = y[i__4].i + 
				q__2.i;
			y[i__3].r = q__1.r, y[i__3].i = q__1.i;
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = jx;
		if (x[i__2].r != (float)0. || x[i__2].i != (float)0.) {
		    i__2 = jx;
		    q__1.r = alpha->r * x[i__2].r - alpha->i * x[i__2].i, 
			    q__1.i = alpha->r * x[i__2].i + alpha->i * x[i__2]
			    .r;
		    temp.r = q__1.r, temp.i = q__1.i;
		    iy = ky;
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = iy;
			i__4 = iy;
			i__5 = i__ + j * a_dim1;
			q__2.r = temp.r * a[i__5].r - temp.i * a[i__5].i, 
				q__2.i = temp.r * a[i__5].i + temp.i * a[i__5]
				.r;
			q__1.r = y[i__4].r + q__2.r, q__1.i = y[i__4].i + 
				q__2.i;
			y[i__3].r = q__1.r, y[i__3].i = q__1.i;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
 */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp.r = (float)0., temp.i = (float)0.;
		if (noconj) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * a_dim1;
			i__4 = i__;
			q__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
				.i, q__2.i = a[i__3].r * x[i__4].i + a[i__3]
				.i * x[i__4].r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L90: */
		    }
		} else {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			r_cnjg(&q__3, &a[i__ + j * a_dim1]);
			i__3 = i__;
			q__2.r = q__3.r * x[i__3].r - q__3.i * x[i__3].i, 
				q__2.i = q__3.r * x[i__3].i + q__3.i * x[i__3]
				.r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
/* L100: */
		    }
		}
		i__2 = jy;
		i__3 = jy;
		q__2.r = alpha->r * temp.r - alpha->i * temp.i, q__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
		q__1.r = y[i__3].r + q__2.r, q__1.i = y[i__3].i + q__2.i;
		y[i__2].r = q__1.r, y[i__2].i = q__1.i;
		jy += *incy;
/* L110: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp.r = (float)0., temp.i = (float)0.;
		ix = kx;
		if (noconj) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			i__3 = i__ + j * a_dim1;
			i__4 = ix;
			q__2.r = a[i__3].r * x[i__4].r - a[i__3].i * x[i__4]
				.i, q__2.i = a[i__3].r * x[i__4].i + a[i__3]
				.i * x[i__4].r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
			ix += *incx;
/* L120: */
		    }
		} else {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			r_cnjg(&q__3, &a[i__ + j * a_dim1]);
			i__3 = ix;
			q__2.r = q__3.r * x[i__3].r - q__3.i * x[i__3].i, 
				q__2.i = q__3.r * x[i__3].i + q__3.i * x[i__3]
				.r;
			q__1.r = temp.r + q__2.r, q__1.i = temp.i + q__2.i;
			temp.r = q__1.r, temp.i = q__1.i;
			ix += *incx;
/* L130: */
		    }
		}
		i__2 = jy;
		i__3 = jy;
		q__2.r = alpha->r * temp.r - alpha->i * temp.i, q__2.i = 
			alpha->r * temp.i + alpha->i * temp.r;
		q__1.r = y[i__3].r + q__2.r, q__1.i = y[i__3].i + q__2.i;
		y[i__2].r = q__1.r, y[i__2].i = q__1.i;
		jy += *incy;
/* L140: */
	    }
	}
    }

    return 0;

/*     End of CGEMV . */

} /* cgemv_ */

