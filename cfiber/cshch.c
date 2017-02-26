/* cshch.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* DECK CSHCH */
/* Subroutine */ int cshch_(complex *z__, complex *csh, complex *cch)
{
    /* System generated locals */
    complex q__1;

    /* Builtin functions */
    double r_imag(), sinh(), cosh(), sin(), cos();

    /* Local variables */
    static real cchi, cchr, cshi, cshr, x, y, ch, cn, sh, sn;

/* ***BEGIN PROLOGUE  CSHCH */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBESH and CBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CSHCH-A, ZSHCH-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     CSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y) */
/*     AND CCH=COSH(X+I*Y), WHERE I**2=-1. */

/* ***SEE ALSO  CBESH, CBESK */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CSHCH */
/* ***FIRST EXECUTABLE STATEMENT  CSHCH */
    x = z__->r;
    y = r_imag(z__);
    sh = sinh(x);
    ch = cosh(x);
    sn = sin(y);
    cn = cos(y);
    cshr = sh * cn;
    cshi = ch * sn;
    q__1.r = cshr, q__1.i = cshi;
    csh->r = q__1.r, csh->i = q__1.i;
    cchr = ch * cn;
    cchi = sh * sn;
    q__1.r = cchr, q__1.i = cchi;
    cch->r = q__1.r, cch->i = q__1.i;
    return 0;
} /* cshch_ */

