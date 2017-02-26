/* cswap.f -- translated by f2c (version 19960717).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* DECK CSWAP */
/* Subroutine */ int cswap_(integer *n, complex *cx, integer *incx,
			    complex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__;
    static complex ctemp;
    static integer ns, kx, ky;

/* ***BEGIN PROLOGUE  CSWAP */
/* ***PURPOSE  Interchange two vectors. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A5 */
/* ***TYPE      COMPLEX (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I) */
/* ***KEYWORDS  BLAS, INTERCHANGE, LINEAR ALGEBRA, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       CX  complex vector with N elements */
/*     INCX  storage spacing between elements of CX */
/*       CY  complex vector with N elements */
/*     INCY  storage spacing between elements of CY */

/*     --Output-- */
/*       CX  input vector CY (unchanged if N .LE. 0) */
/*       CY  input vector CX (unchanged if N .LE. 0) */

/*     Interchange complex CX and complex CY */
/*     For I = 0 to N-1, interchange  CX(LX+I*INCX) and CY(LY+I*INCY), */
/*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/*     defined in a similar way using INCY. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical 
*/
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CSWAP */
/* ***FIRST EXECUTABLE STATEMENT  CSWAP */
    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == *incy && *incx > 0) {
	goto L20;
    }

/*     Code for unequal or nonpositive increments. */

    kx = 1;
    ky = 1;
    if (*incx < 0) {
	kx = (1 - *n) * *incx + 1;
    }
    if (*incy < 0) {
	ky = (1 - *n) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = kx;
	ctemp.r = cx[i__2].r, ctemp.i = cx[i__2].i;
	i__2 = kx;
	i__3 = ky;
	cx[i__2].r = cy[i__3].r, cx[i__2].i = cy[i__3].i;
	i__2 = ky;
	cy[i__2].r = ctemp.r, cy[i__2].i = ctemp.i;
	kx += *incx;
	ky += *incy;
/* L10: */
    }
    return 0;

/*     Code for equal, positive, non-unit increments. */

L20:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	i__3 = i__;
	ctemp.r = cx[i__3].r, ctemp.i = cx[i__3].i;
	i__3 = i__;
	i__4 = i__;
	cx[i__3].r = cy[i__4].r, cx[i__3].i = cy[i__4].i;
	i__3 = i__;
	cy[i__3].r = ctemp.r, cy[i__3].i = ctemp.i;
/* L30: */
    }
    return 0;
} /* cswap_ */

