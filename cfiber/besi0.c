/* besi0.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__12 = 12;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK BESI0 */
doublereal besi0_(x)
real *x;
{
    /* Initialized data */

    static real bi0cs[12] = { (float)-.07660547252839144951,(float)
	    1.92733795399380827,(float).2282644586920301339,(float)
	    .01304891466707290428,(float)4.3442709008164874e-4,(float)
	    9.42265768600193e-6,(float)1.4340062895106e-7,(float)
	    1.61384906966e-9,(float)1.396650044e-11,(float)9.579451e-14,(
	    float)5.3339e-16,(float)2.45e-18 };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(), log(), exp();

    /* Local variables */
    static real xmax, xsml, y;
    extern doublereal csevl_();
    extern integer inits_();
    extern doublereal besi0e_(), r1mach_();
    extern /* Subroutine */ int xermsg_();
    static integer nti0;

/* ***BEGIN PROLOGUE  BESI0 */
/* ***PURPOSE  Compute the hyperbolic Bessel function of the first kind */
/*            of order zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESI0-S, DBESI0-D) */
/* ***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESI0(X) computes the modified (hyperbolic) Bessel function */
/* of the first kind of order zero and real argument X. */

/* Series for BI0        on the interval  0.          to  9.00000D+00 */
/*                                        with weighted error   2.46E-18 
*/
/*                                         log weighted error  17.61 */
/*                               significant figures required  17.90 */
/*                                    decimal places required  18.15 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI0E, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESI0 */
/* ***FIRST EXECUTABLE STATEMENT  BESI0 */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	nti0 = inits_(bi0cs, &c__12, &r__1);
	xsml = sqrt(r1mach_(&c__3) * (float)4.5);
	xmax = log(r1mach_(&c__2));
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > (float)3.) {
	goto L20;
    }

    ret_val = (float)1.;
    if (y > xsml) {
	r__1 = y * y / (float)4.5 - (float)1.;
	ret_val = csevl_(&r__1, bi0cs, &nti0) + (float)2.75;
    }
    return ret_val;

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "BESI0", "ABS(X) SO BIG I0 OVERFLOWS", &c__1, &c__2,
		 6L, 5L, 26L);
    }

    ret_val = exp(y) * besi0e_(x);

    return ret_val;
} /* besi0_ */

