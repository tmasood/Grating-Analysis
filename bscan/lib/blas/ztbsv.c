
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int ztbsv_(char *uplo, char *trans, char *diag, integer *n, 
	integer *k, doublecomplex *a, integer *lda, doublecomplex *x, integer 
	*incx)
{


    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer info;
    static doublecomplex temp;
    static integer i, j, l;
    extern logical lsame_(char *, char *);
    static integer kplus1, ix, jx, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical noconj, nounit;


/*  Purpose   
    =======   

    ZTBSV  solves one of the systems of equations   

       A*x = b,   or   A'*x = b,   or   conjg( A' )*x = b,   

    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular band matrix, with ( k + 1 )   
    diagonals.   

    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   

                TRANS = 'N' or 'n'   A*x = b.   

                TRANS = 'T' or 't'   A'*x = b.   

                TRANS = 'C' or 'c'   conjg( A' )*x = b.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    K      - INTEGER.   
             On entry with UPLO = 'U' or 'u', K specifies the number of   
             super-diagonals of the matrix A.   
             On entry with UPLO = 'L' or 'l', K specifies the number of   
             sub-diagonals of the matrix A.   
             K must satisfy  0 .le. K.   
             Unchanged on exit.   

    A      - COMPLEX*16       array of DIMENSION ( LDA, n ).   
             Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )   
             by n part of the array A must contain the upper triangular   
             band part of the matrix of coefficients, supplied column by 
  
             column, with the leading diagonal of the matrix in row   
             ( k + 1 ) of the array, the first super-diagonal starting at 
  
             position 2 in row k, and so on. The top left k by k triangle 
  
             of the array A is not referenced.   
             The following program segment will transfer an upper   
             triangular band matrix from conventional full matrix storage 
  
             to band storage:   

                   DO 20, J = 1, N   
                      M = K + 1 - J   
                      DO 10, I = MAX( 1, J - K ), J   
                         A( M + I, J ) = matrix( I, J )   
                10    CONTINUE   
                20 CONTINUE   

             Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )   
             by n part of the array A must contain the lower triangular   
             band part of the matrix of coefficients, supplied column by 
  
             column, with the leading diagonal of the matrix in row 1 of 
  
             the array, the first sub-diagonal starting at position 1 in 
  
             row 2, and so on. The bottom right k by k triangle of the   
             array A is not referenced.   
             The following program segment will transfer a lower   
             triangular band matrix from conventional full matrix storage 
  
             to band storage:   

                   DO 20, J = 1, N   
                      M = 1 - J   
                      DO 10, I = J, MIN( N, J + K )   
                         A( M + I, J ) = matrix( I, J )   
                10    CONTINUE   
                20 CONTINUE   

             Note that when DIAG = 'U' or 'u' the elements of the array A 
  
             corresponding to the diagonal elements of the matrix are not 
  
             referenced, but are assumed to be unity.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             ( k + 1 ).   
             Unchanged on exit.   

    X      - COMPLEX*16       array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten 
  
             with the solution vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
	info = 2;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < *k + 1) {
	info = 7;
    } else if (*incx == 0) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("ZTBSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    noconj = lsame_(trans, "T");
    nounit = lsame_(diag, "N");

/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are   
       accessed by sequentially with one pass through A. */

    if (lsame_(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U")) {
	    kplus1 = *k + 1;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    i__1 = j;
		    if (X(j).r != 0. || X(j).i != 0.) {
			l = kplus1 - j;
			if (nounit) {
			    i__1 = j;
			    z_div(&z__1, &X(j), &A(kplus1,j));
			    X(j).r = z__1.r, X(j).i = z__1.i;
			}
			i__1 = j;
			temp.r = X(j).r, temp.i = X(j).i;
/* Computing MAX */
			i__2 = 1, i__3 = j - *k;
			i__1 = max(i__2,i__3);
			for (i = j - 1; i >= max(1,j-*k); --i) {
			    i__2 = i;
			    i__3 = i;
			    i__4 = l + i + j * a_dim1;
			    z__2.r = temp.r * A(l+i,j).r - temp.i * A(l+i,j).i, 
				    z__2.i = temp.r * A(l+i,j).i + temp.i * A(l+i,j).r;
			    z__1.r = X(i).r - z__2.r, z__1.i = X(i).i - 
				    z__2.i;
			    X(i).r = z__1.r, X(i).i = z__1.i;
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    kx -= *incx;
		    i__1 = jx;
		    if (X(jx).r != 0. || X(jx).i != 0.) {
			ix = kx;
			l = kplus1 - j;
			if (nounit) {
			    i__1 = jx;
			    z_div(&z__1, &X(jx), &A(kplus1,j));
			    X(jx).r = z__1.r, X(jx).i = z__1.i;
			}
			i__1 = jx;
			temp.r = X(jx).r, temp.i = X(jx).i;
/* Computing MAX */
			i__2 = 1, i__3 = j - *k;
			i__1 = max(i__2,i__3);
			for (i = j - 1; i >= max(1,j-*k); --i) {
			    i__2 = ix;
			    i__3 = ix;
			    i__4 = l + i + j * a_dim1;
			    z__2.r = temp.r * A(l+i,j).r - temp.i * A(l+i,j).i, 
				    z__2.i = temp.r * A(l+i,j).i + temp.i * A(l+i,j).r;
			    z__1.r = X(ix).r - z__2.r, z__1.i = X(ix).i - 
				    z__2.i;
			    X(ix).r = z__1.r, X(ix).i = z__1.i;
			    ix -= *incx;
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    if (X(j).r != 0. || X(j).i != 0.) {
			l = 1 - j;
			if (nounit) {
			    i__2 = j;
			    z_div(&z__1, &X(j), &A(1,j));
			    X(j).r = z__1.r, X(j).i = z__1.i;
			}
			i__2 = j;
			temp.r = X(j).r, temp.i = X(j).i;
/* Computing MIN */
			i__3 = *n, i__4 = j + *k;
			i__2 = min(i__3,i__4);
			for (i = j + 1; i <= min(*n,j+*k); ++i) {
			    i__3 = i;
			    i__4 = i;
			    i__5 = l + i + j * a_dim1;
			    z__2.r = temp.r * A(l+i,j).r - temp.i * A(l+i,j).i, 
				    z__2.i = temp.r * A(l+i,j).i + temp.i * A(l+i,j).r;
			    z__1.r = X(i).r - z__2.r, z__1.i = X(i).i - 
				    z__2.i;
			    X(i).r = z__1.r, X(i).i = z__1.i;
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    kx += *incx;
		    i__2 = jx;
		    if (X(jx).r != 0. || X(jx).i != 0.) {
			ix = kx;
			l = 1 - j;
			if (nounit) {
			    i__2 = jx;
			    z_div(&z__1, &X(jx), &A(1,j));
			    X(jx).r = z__1.r, X(jx).i = z__1.i;
			}
			i__2 = jx;
			temp.r = X(jx).r, temp.i = X(jx).i;
/* Computing MIN */
			i__3 = *n, i__4 = j + *k;
			i__2 = min(i__3,i__4);
			for (i = j + 1; i <= min(*n,j+*k); ++i) {
			    i__3 = ix;
			    i__4 = ix;
			    i__5 = l + i + j * a_dim1;
			    z__2.r = temp.r * A(l+i,j).r - temp.i * A(l+i,j).i, 
				    z__2.i = temp.r * A(l+i,j).i + temp.i * A(l+i,j).r;
			    z__1.r = X(ix).r - z__2.r, z__1.i = X(ix).i - 
				    z__2.i;
			    X(ix).r = z__1.r, X(ix).i = z__1.i;
			    ix += *incx;
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x  or  x := inv( conjg( A') )*x. */

	if (lsame_(uplo, "U")) {
	    kplus1 = *k + 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__2 = j;
		    temp.r = X(j).r, temp.i = X(j).i;
		    l = kplus1 - j;
		    if (noconj) {
/* Computing MAX */
			i__2 = 1, i__3 = j - *k;
			i__4 = j - 1;
			for (i = max(1,j-*k); i <= j-1; ++i) {
			    i__2 = l + i + j * a_dim1;
			    i__3 = i;
			    z__2.r = A(l+i,j).r * X(i).r - A(l+i,j).i * X(
				    i).i, z__2.i = A(l+i,j).r * X(i).i + 
				    A(l+i,j).i * X(i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
/* L90: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &A(kplus1,j));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
/* Computing MAX */
			i__4 = 1, i__2 = j - *k;
			i__3 = j - 1;
			for (i = max(1,j-*k); i <= j-1; ++i) {
			    d_cnjg(&z__3, &A(l+i,j));
			    i__4 = i;
			    z__2.r = z__3.r * X(i).r - z__3.i * X(i).i, 
				    z__2.i = z__3.r * X(i).i + z__3.i * X(
				    i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
/* L100: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &A(kplus1,j));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__3 = j;
		    X(j).r = temp.r, X(j).i = temp.i;
/* L110: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    i__3 = jx;
		    temp.r = X(jx).r, temp.i = X(jx).i;
		    ix = kx;
		    l = kplus1 - j;
		    if (noconj) {
/* Computing MAX */
			i__3 = 1, i__4 = j - *k;
			i__2 = j - 1;
			for (i = max(1,j-*k); i <= j-1; ++i) {
			    i__3 = l + i + j * a_dim1;
			    i__4 = ix;
			    z__2.r = A(l+i,j).r * X(ix).r - A(l+i,j).i * X(
				    ix).i, z__2.i = A(l+i,j).r * X(ix).i + 
				    A(l+i,j).i * X(ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix += *incx;
/* L120: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &A(kplus1,j));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
/* Computing MAX */
			i__2 = 1, i__3 = j - *k;
			i__4 = j - 1;
			for (i = max(1,j-*k); i <= j-1; ++i) {
			    d_cnjg(&z__3, &A(l+i,j));
			    i__2 = ix;
			    z__2.r = z__3.r * X(ix).r - z__3.i * X(ix).i, 
				    z__2.i = z__3.r * X(ix).i + z__3.i * X(
				    ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix += *incx;
/* L130: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &A(kplus1,j));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__4 = jx;
		    X(jx).r = temp.r, X(jx).i = temp.i;
		    jx += *incx;
		    if (j > *k) {
			kx += *incx;
		    }
/* L140: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    i__1 = j;
		    temp.r = X(j).r, temp.i = X(j).i;
		    l = 1 - j;
		    if (noconj) {
/* Computing MIN */
			i__1 = *n, i__4 = j + *k;
			i__2 = j + 1;
			for (i = min(*n,j+*k); i >= j+1; --i) {
			    i__1 = l + i + j * a_dim1;
			    i__4 = i;
			    z__2.r = A(l+i,j).r * X(i).r - A(l+i,j).i * X(
				    i).i, z__2.i = A(l+i,j).r * X(i).i + 
				    A(l+i,j).i * X(i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
/* L150: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &A(1,j));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
/* Computing MIN */
			i__2 = *n, i__1 = j + *k;
			i__4 = j + 1;
			for (i = min(*n,j+*k); i >= j+1; --i) {
			    d_cnjg(&z__3, &A(l+i,j));
			    i__2 = i;
			    z__2.r = z__3.r * X(i).r - z__3.i * X(i).i, 
				    z__2.i = z__3.r * X(i).i + z__3.i * X(
				    i).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
/* L160: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &A(1,j));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__4 = j;
		    X(j).r = temp.r, X(j).i = temp.i;
/* L170: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    i__4 = jx;
		    temp.r = X(jx).r, temp.i = X(jx).i;
		    ix = kx;
		    l = 1 - j;
		    if (noconj) {
/* Computing MIN */
			i__4 = *n, i__2 = j + *k;
			i__1 = j + 1;
			for (i = min(*n,j+*k); i >= j+1; --i) {
			    i__4 = l + i + j * a_dim1;
			    i__2 = ix;
			    z__2.r = A(l+i,j).r * X(ix).r - A(l+i,j).i * X(
				    ix).i, z__2.i = A(l+i,j).r * X(ix).i + 
				    A(l+i,j).i * X(ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix -= *incx;
/* L180: */
			}
			if (nounit) {
			    z_div(&z__1, &temp, &A(1,j));
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    } else {
/* Computing MIN */
			i__1 = *n, i__4 = j + *k;
			i__2 = j + 1;
			for (i = min(*n,j+*k); i >= j+1; --i) {
			    d_cnjg(&z__3, &A(l+i,j));
			    i__1 = ix;
			    z__2.r = z__3.r * X(ix).r - z__3.i * X(ix).i, 
				    z__2.i = z__3.r * X(ix).i + z__3.i * X(
				    ix).r;
			    z__1.r = temp.r - z__2.r, z__1.i = temp.i - 
				    z__2.i;
			    temp.r = z__1.r, temp.i = z__1.i;
			    ix -= *incx;
/* L190: */
			}
			if (nounit) {
			    d_cnjg(&z__2, &A(1,j));
			    z_div(&z__1, &temp, &z__2);
			    temp.r = z__1.r, temp.i = z__1.i;
			}
		    }
		    i__2 = jx;
		    X(jx).r = temp.r, X(jx).i = temp.i;
		    jx -= *incx;
		    if (*n - j >= *k) {
			kx -= *incx;
		    }
/* L200: */
		}
	    }
	}
    }

    return 0;

/*     End of ZTBSV . */

} /* ztbsv_ */

