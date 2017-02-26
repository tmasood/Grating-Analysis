#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static integer c__4 = 4;
static integer c__3 = 3;
static integer c__1 = 1;

/* DECK ALNGAM */
doublereal alngam_(real *x)
{
    /* Initialized data */

    static real sq2pil = (float).91893853320467274;
    static real sqpi2l = (float).22579135264472743;
    static real pi = (float)3.14159265358979324;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double log(), sqrt(), sin(), r_int();

    /* Local variables */
    static real xmax;
    extern doublereal gamma_();
    static real y, dxrel;
    extern doublereal r1mach_(), r9lgmc_();
    extern /* Subroutine */ int xermsg_();
    static real sinpiy;

/* ***BEGIN PROLOGUE  ALNGAM */
/* ***PURPOSE  Compute the logarithm of the absolute value of the Gamma */
/*            function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      SINGLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C) */
/* ***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
 */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* ALNGAM(X) computes the logarithm of the absolute value of the */
/* gamma function at X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  GAMMA, R1MACH, R9LGMC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  ALNGAM */
/* ***FIRST EXECUTABLE STATEMENT  ALNGAM */
    if (first) {
	xmax = r1mach_(&c__2) / log(r1mach_(&c__2));
	dxrel = sqrt(r1mach_(&c__4));
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > (float)10.) {
	goto L20;
    }

/* LOG (ABS (GAMMA(X))) FOR  ABS(X) .LE. 10.0 */

    ret_val = log((r__1 = gamma_(x), dabs(r__1)));

    return ret_val;

/* LOG (ABS (GAMMA(X))) FOR ABS(X) .GT. 10.0 */

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "ALNGAM", "ABS(X) SO BIG ALNGAM OVERFLOWS", &c__2, &
		c__2, 6L, 6L, 30L);
    }

    if (*x > (float)0.) {
	ret_val = sq2pil + (*x - (float).5) * log(*x) - *x + r9lgmc_(&y);
    }
    if (*x > (float)0.) {
	return ret_val;
    }

    sinpiy = (r__1 = sin(pi * y), dabs(r__1));
    if (sinpiy == (float)0.) {
	xermsg_("SLATEC", "ALNGAM", "X IS A NEGATIVE INTEGER", &c__3, &c__2, 
		6L, 6L, 23L);
    }

    r__2 = *x - (float).5;
    if ((r__1 = (*x - r_int(&r__2)) / *x, dabs(r__1)) < dxrel) {
	xermsg_("SLATEC", "ALNGAM", "ANSWER LT HALF PRECISION BECAUSE X TOO \
NEAR NEGATIVE INTEGER", &c__1, &c__1, 6L, 6L, 60L);
    }

    ret_val = sqpi2l + (*x - (float).5) * log(y) - *x - log(sinpiy) - r9lgmc_(
	    &y);
    return ret_val;

} /* alngam_ */

