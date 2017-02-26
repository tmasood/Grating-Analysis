/* besi1.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__11 = 11;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK BESI1 */
doublereal besi1_(x)
real *x;
{
    /* Initialized data */

    static real bi1cs[11] = { (float)-.001971713261099859,(float)
	    .40734887667546481,(float).034838994299959456,(float)
	    .001545394556300123,(float)4.1888521098377e-5,(float)
	    7.64902676483e-7,(float)1.0042493924e-8,(float)9.9322077e-11,(
	    float)7.6638e-13,(float)4.741e-15,(float)2.4e-17 };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(), log(), exp();

    /* Local variables */
    static real xmin, xmax, xsml, y;
    extern doublereal csevl_();
    extern integer inits_();
    extern doublereal besi1e_(), r1mach_();
    extern /* Subroutine */ int xermsg_();
    static integer nti1;

/* ***BEGIN PROLOGUE  BESI1 */
/* ***PURPOSE  Compute the modified (hyperbolic) Bessel function of the */
/*            first kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESI1-S, DBESI1-D) */
/* ***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESI1(X) calculates the modified (hyperbolic) Bessel function */
/* of the first kind of order one for real argument X. */

/* Series for BI1        on the interval  0.          to  9.00000D+00 */
/*                                        with weighted error   2.40E-17 
*/
/*                                         log weighted error  16.62 */
/*                               significant figures required  16.23 */
/*                                    decimal places required  17.14 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI1E, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESI1 */
/* ***FIRST EXECUTABLE STATEMENT  BESI1 */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	nti1 = inits_(bi1cs, &c__11, &r__1);
	xmin = r1mach_(&c__1) * (float)2.;
	xsml = sqrt(r1mach_(&c__3) * (float)4.5);
	xmax = log(r1mach_(&c__2));
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > (float)3.) {
	goto L20;
    }

    ret_val = (float)0.;
    if (y == (float)0.) {
	return ret_val;
    }

    if (y <= xmin) {
	xermsg_("SLATEC", "BESI1", "ABS(X) SO SMALL I1 UNDERFLOWS", &c__1, &
		c__1, 6L, 5L, 29L);
    }
    if (y > xmin) {
	ret_val = *x * (float).5;
    }
    if (y > xsml) {
	r__1 = y * y / (float)4.5 - (float)1.;
	ret_val = *x * (csevl_(&r__1, bi1cs, &nti1) + (float).875);
    }
    return ret_val;

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "BESI1", "ABS(X) SO BIG I1 OVERFLOWS", &c__2, &c__2,
		 6L, 5L, 26L);
    }

    ret_val = exp(y) * besi1e_(x);

    return ret_val;
} /* besi1_ */

