/* besk0e.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__11 = 11;
static integer c__3 = 3;
static integer c__17 = 17;
static integer c__14 = 14;
static integer c__2 = 2;

/* DECK BESK0E */
doublereal besk0e_(x)
real *x;
{
    /* Initialized data */

    static real bk0cs[11] = { (float)-.03532739323390276872,(float)
	    .3442898999246284869,(float).03597993651536150163,(float)
	    .00126461541144692592,(float)2.286212103119451e-5,(float)
	    2.5347910790261e-7,(float)1.90451637722e-9,(float)1.034969525e-11,
	    (float)4.259816e-14,(float)1.3744e-16,(float)3.5e-19 };
    static real ak0cs[17] = { (float)-.07643947903327941,(float)
	    -.02235652605699819,(float)7.7341811546938e-4,(float)
	    -4.281006688886e-5,(float)3.08170017386e-6,(float)-2.639367222e-7,
	    (float)2.563713036e-8,(float)-2.74270554e-9,(float)3.1694296e-10,(
	    float)-3.902353e-11,(float)5.06804e-12,(float)-6.8895e-13,(float)
	    9.744e-14,(float)-1.427e-14,(float)2.15e-15,(float)-3.3e-16,(
	    float)5e-17 };
    static real ak02cs[14] = { (float)-.01201869826307592,(float)
	    -.00917485269102569,(float)1.444550931775e-4,(float)
	    -4.01361417543e-6,(float)1.5678318108e-7,(float)-7.77011043e-9,(
	    float)4.6111825e-10,(float)-3.158592e-11,(float)2.43501e-12,(
	    float)-2.0743e-13,(float)1.925e-14,(float)-1.92e-15,(float)2e-16,(
	    float)-2e-17 };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(), exp(), log();

    /* Local variables */
    static real xsml;
    extern doublereal besi0_();
    static integer ntak0;
    static real y;
    static integer ntak02;
    extern doublereal csevl_();
    extern integer inits_();
    extern doublereal r1mach_();
    extern /* Subroutine */ int xermsg_();
    static integer ntk0;

/* ***BEGIN PROLOGUE  BESK0E */
/* ***PURPOSE  Compute the exponentially scaled modified (hyperbolic) */
/*            Bessel function of the third kind of order zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESK0E-S, DBSK0E-D) */
/* ***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION, 
*/
/*             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESK0E(X) computes the exponentially scaled modified (hyperbolic) */
/* Bessel function of third kind of order zero for real argument */
/* X .GT. 0.0, i.e., EXP(X)*K0(X). */

/* Series for BK0        on the interval  0.          to  4.00000D+00 */
/*                                        with weighted error   3.57E-19 
*/
/*                                         log weighted error  18.45 */
/*                               significant figures required  17.99 */
/*                                    decimal places required  18.97 */

/* Series for AK0        on the interval  1.25000D-01 to  5.00000D-01 */
/*                                        with weighted error   5.34E-17 
*/
/*                                         log weighted error  16.27 */
/*                               significant figures required  14.92 */
/*                                    decimal places required  16.89 */

/* Series for AK02       on the interval  0.          to  1.25000D-01 */
/*                                        with weighted error   2.34E-17 
*/
/*                                         log weighted error  16.63 */
/*                               significant figures required  14.67 */
/*                                    decimal places required  17.20 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI0, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESK0E */
/* ***FIRST EXECUTABLE STATEMENT  BESK0E */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	ntk0 = inits_(bk0cs, &c__11, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntak0 = inits_(ak0cs, &c__17, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntak02 = inits_(ak02cs, &c__14, &r__1);
	xsml = sqrt(r1mach_(&c__3) * (float)4.);
    }
    first = FALSE_;

    if (*x <= (float)0.) {
	xermsg_("SLATEC", "BESK0E", "X IS ZERO OR NEGATIVE", &c__2, &c__2, 6L,
		 6L, 21L);
    }
    if (*x > (float)2.) {
	goto L20;
    }

    y = (float)0.;
    if (*x > xsml) {
	y = *x * *x;
    }
    r__1 = y * (float).5 - (float)1.;
    ret_val = exp(*x) * (-log(*x * (float).5) * besi0_(x) - (float).25 + 
	    csevl_(&r__1, bk0cs, &ntk0));
    return ret_val;

L20:
    if (*x <= (float)8.) {
	r__1 = ((float)16. / *x - (float)5.) / (float)3.;
	ret_val = (csevl_(&r__1, ak0cs, &ntak0) + (float)1.25) / sqrt(*x);
    }
    if (*x > (float)8.) {
	r__1 = (float)16. / *x - (float)1.;
	ret_val = (csevl_(&r__1, ak02cs, &ntak02) + (float)1.25) / sqrt(*x);
    }

    return ret_val;
} /* besk0e_ */

