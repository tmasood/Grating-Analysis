/* besk1.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__11 = 11;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK BESK1 */
doublereal besk1_(x)
real *x;
{
    /* Initialized data */

    static real bk1cs[11] = { (float).0253002273389477705,(float)
	    -.353155960776544876,(float)-.122611180822657148,(float)
	    -.0069757238596398643,(float)-1.730288957513052e-4,(float)
	    -2.4334061415659e-6,(float)-2.21338763073e-8,(float)
	    -1.411488392e-10,(float)-6.666901e-13,(float)-2.4274e-15,(float)
	    -7e-18 };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double log(), exp(), sqrt();

    /* Local variables */
    static real xmin, xmax, xsml;
    extern doublereal besi1_();
    static real y;
    extern doublereal csevl_();
    extern integer inits_();
    static real xmaxt;
    extern doublereal besk1e_(), r1mach_();
    extern /* Subroutine */ int xermsg_();
    static integer ntk1;

/* ***BEGIN PROLOGUE  BESK1 */
/* ***PURPOSE  Compute the modified (hyperbolic) Bessel function of the */
/*            third kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESK1-S, DBESK1-D) */
/* ***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESK1(X) computes the modified (hyperbolic) Bessel function of third */
/* kind of order one for real argument X, where X .GT. 0. */

/* Series for BK1        on the interval  0.          to  4.00000D+00 */
/*                                        with weighted error   7.02E-18 
*/
/*                                         log weighted error  17.15 */
/*                               significant figures required  16.73 */
/*                                    decimal places required  17.67 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI1, BESK1E, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESK1 */
/* ***FIRST EXECUTABLE STATEMENT  BESK1 */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	ntk1 = inits_(bk1cs, &c__11, &r__1);
/* Computing MAX */
	r__1 = log(r1mach_(&c__1)), r__2 = -log(r1mach_(&c__2));
	xmin = exp(dmax(r__1,r__2) + (float).01);
	xsml = sqrt(r1mach_(&c__3) * (float)4.);
	xmaxt = -log(r1mach_(&c__1));
	xmax = xmaxt - xmaxt * (float).5 * log(xmaxt) / (xmaxt + (float).5);
    }
    first = FALSE_;

    if (*x <= (float)0.) {
	xermsg_("SLATEC", "BESK1", "X IS ZERO OR NEGATIVE", &c__2, &c__2, 6L, 
		5L, 21L);
    }
    if (*x > (float)2.) {
	goto L20;
    }

    if (*x < xmin) {
	xermsg_("SLATEC", "BESK1", "X SO SMALL K1 OVERFLOWS", &c__3, &c__2, 
		6L, 5L, 23L);
    }
    y = (float)0.;
    if (*x > xsml) {
	y = *x * *x;
    }
    r__1 = y * (float).5 - (float)1.;
    ret_val = log(*x * (float).5) * besi1_(x) + (csevl_(&r__1, bk1cs, &ntk1) 
	    + (float).75) / *x;
    return ret_val;

L20:
    ret_val = (float)0.;
    if (*x > xmax) {
	xermsg_("SLATEC", "BESK1", "X SO BIG K1 UNDERFLOWS", &c__1, &c__1, 6L,
		 5L, 22L);
    }
    if (*x > xmax) {
	return ret_val;
    }

    ret_val = exp(-(*x)) * besk1e_(x);

    return ret_val;
} /* besk1_ */

