/* cbunk.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* DECK CBUNK */
/* Subroutine */ int cbunk_(complex *z__, real *fnu, integer *kode,
			    integer *mr, integer *n, complex *y,
			    integer *nz, real *tol, real *elim,
			    real *alim)
{
    /* Builtin functions */
    double r_imag();

    /* Local variables */
    extern /* Subroutine */ int cunk1_(complex *, real *, integer *,
			    integer *, integer *, complex *,
			    integer *, real *, real *, real *);
    extern int cunk2_(complex *, real *, integer *, integer *, integer *,
		      complex *, integer *, real *, real *, real *);
    static real ax, ay, xx, yy;

/* ***BEGIN PROLOGUE  CBUNK */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBESH and CBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CBUNK-A, ZBUNK-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     CBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU.GT.FNUL. */
/*     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z) */
/*     IN CUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN CUNK2 */

/* ***SEE ALSO  CBESH, CBESK */
/* ***ROUTINES CALLED  CUNK1, CUNK2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CBUNK */
/* ***FIRST EXECUTABLE STATEMENT  CBUNK */
    /* Parameter adjustments */
    --y;

    /* Function Body */
    *nz = 0;
    xx = z__->r;
    yy = r_imag(z__);
    ax = dabs(xx) * (float)1.7321;
    ay = dabs(yy);
    if (ay > ax) {
	goto L10;
    }
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN */
/*     -PI/3.LE.ARG(Z).LE.PI/3 */
/* ----------------------------------------------------------------------- */
    cunk1_(z__, fnu, kode, mr, n, &y[1], nz, tol, elim, alim);
    goto L20;
L10:
/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU */
/*     APPLIED IN PI/3.LT.ABS(ARG(Z)).LE.PI/2 WHERE M=+I OR -I */
/*     AND HPI=PI/2 */
/* ----------------------------------------------------------------------- */
    cunk2_(z__, fnu, kode, mr, n, &y[1], nz, tol, elim, alim);
L20:
    return 0;
} /* cbunk_ */

