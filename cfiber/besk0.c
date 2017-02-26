/* besk0.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__11 = 11;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK BESK0 */
doublereal besk0_(x)
real *x;
{
    /* Initialized data */

    static real bk0cs[11] = { (float)-.03532739323390276872,(float)
	    .3442898999246284869,(float).03597993651536150163,(float)
	    .00126461541144692592,(float)2.286212103119451e-5,(float)
	    2.5347910790261e-7,(float)1.90451637722e-9,(float)1.034969525e-11,
	    (float)4.259816e-14,(float)1.3744e-16,(float)3.5e-19 };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(), log(), exp();

    /* Local variables */
    static real xmax, xsml;
    extern doublereal besi0_();
    static real y;
    extern doublereal csevl_();
    extern integer inits_();
    static real xmaxt;
    extern doublereal besk0e_(), r1mach_();
    extern /* Subroutine */ int xermsg_();
    static integer ntk0;

/* ***BEGIN PROLOGUE  BESK0 */
/* ***PURPOSE  Compute the modified (hyperbolic) Bessel function of the */
/*            third kind of order zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESK0-S, DBESK0-D) */
/* ***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESK0(X) calculates the modified (hyperbolic) Bessel function */
/* of the third kind of order zero for real argument X .GT. 0.0. */

/* Series for BK0        on the interval  0.          to  4.00000D+00 */
/*                                        with weighted error   3.57E-19 
*/
/*                                         log weighted error  18.45 */
/*                               significant figures required  17.99 */
/*                                    decimal places required  18.97 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI0, BESK0E, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESK0 */
/* ***FIRST EXECUTABLE STATEMENT  BESK0 */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	ntk0 = inits_(bk0cs, &c__11, &r__1);
	xsml = sqrt(r1mach_(&c__3) * (float)4.);
	xmaxt = -log(r1mach_(&c__1));
	xmax = xmaxt - xmaxt * (float).5 * log(xmaxt) / (xmaxt + (float).5) - 
		(float).01;
    }
    first = FALSE_;

    if (*x <= (float)0.) {
	xermsg_("SLATEC", "BESK0", "X IS ZERO OR NEGATIVE", &c__2, &c__2, 6L, 
		5L, 21L);
    }
    if (*x > (float)2.) {
	goto L20;
    }

    y = (float)0.;
    if (*x > xsml) {
	y = *x * *x;
    }
    r__1 = y * (float).5 - (float)1.;
    ret_val = -log(*x * (float).5) * besi0_(x) - (float).25 + csevl_(&r__1, 
	    bk0cs, &ntk0);
    return ret_val;

L20:
    ret_val = (float)0.;
    if (*x > xmax) {
	xermsg_("SLATEC", "BESK0", "X SO BIG K0 UNDERFLOWS", &c__1, &c__1, 6L,
		 5L, 22L);
    }
    if (*x > xmax) {
	return ret_val;
    }

    ret_val = exp(-(*x)) * besk0e_(x);

    return ret_val;
} /* besk0_ */

