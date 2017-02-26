/* inits.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* DECK INITS */
integer inits_(os, nos, eta)
real *os;
integer *nos;
real *eta;
{
    /* System generated locals */
    integer ret_val, i__1;
    real r__1;

    /* Local variables */
    static integer i__, ii;
    extern /* Subroutine */ int xermsg_();
    static real err;

/* ***BEGIN PROLOGUE  INITS */
/* ***PURPOSE  Determine the number of terms needed in an orthogonal */
/*            polynomial series so that it meets a specified accuracy. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C3A2 */
/* ***TYPE      SINGLE PRECISION (INITS-S, INITDS-D) */
/* ***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL, */
/*             ORTHOGONAL SERIES, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/*  Initialize the orthogonal series, represented by the array OS, so */
/*  that INITS is the number of terms needed to insure the error is no */
/*  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth */
/*  machine precision. */

/*             Input Arguments -- */
/*   OS     single precision array of NOS coefficients in an orthogonal */
/*          series. */
/*   NOS    number of coefficients in OS. */
/*   ETA    single precision scalar containing requested accuracy of */
/*          series. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891115  Modified error message.  (WRB) */
/*   891115  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  INITS */
/* ***FIRST EXECUTABLE STATEMENT  INITS */
    /* Parameter adjustments */
    --os;

    /* Function Body */
    if (*nos < 1) {
	xermsg_("SLATEC", "INITS", "Number of coefficients is less than 1", &
		c__2, &c__1, 6L, 5L, 37L);
    }

    err = (float)0.;
    i__1 = *nos;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *nos + 1 - ii;
	err += (r__1 = os[i__], dabs(r__1));
	if (err > *eta) {
	    goto L20;
	}
/* L10: */
    }

L20:
    if (i__ == *nos) {
	xermsg_("SLATEC", "INITS", "Chebyshev series too short for specified\
 accuracy", &c__1, &c__1, 6L, 5L, 49L);
    }
    ret_val = i__;

    return ret_val;
} /* inits_ */

