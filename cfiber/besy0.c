/* besy0.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__13 = 13;
static integer c__3 = 3;
static integer c__21 = 21;
static integer c__24 = 24;
static integer c__4 = 4;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK BESY0 */
doublereal besy0_(x)
real *x;
{
    /* Initialized data */

    static real by0cs[13] = { (float)-.011277839392865573,(float)
	    -.12834523756042035,(float)-.10437884799794249,(float)
	    .023662749183969695,(float)-.002090391647700486,(float)
	    1.03975453939057e-4,(float)-3.369747162423e-6,(float)
	    7.7293842676e-8,(float)-1.324976772e-9,(float)1.7648232e-11,(
	    float)-1.88105e-13,(float)1.641e-15,(float)-1.1e-17 };
    static real bm0cs[21] = { (float).09284961637381644,(float)
	    -.00142987707403484,(float)2.830579271257e-5,(float)
	    -1.43300611424e-6,(float)1.2028628046e-7,(float)-1.397113013e-8,(
	    float)2.04076188e-9,(float)-3.5399669e-10,(float)7.024759e-11,(
	    float)-1.554107e-11,(float)3.76226e-12,(float)-9.8282e-13,(float)
	    2.7408e-13,(float)-8.091e-14,(float)2.511e-14,(float)-8.14e-15,(
	    float)2.75e-15,(float)-9.6e-16,(float)3.4e-16,(float)-1.2e-16,(
	    float)4e-17 };
    static real bth0cs[24] = { (float)-.24639163774300119,(float)
	    .001737098307508963,(float)-6.2183633402968e-5,(float)
	    4.368050165742e-6,(float)-4.56093019869e-7,(float)6.2197400101e-8,
	    (float)-1.0300442889e-8,(float)1.979526776e-9,(float)
	    -4.28198396e-10,(float)1.0203584e-10,(float)-2.6363898e-11,(float)
	    7.297935e-12,(float)-2.144188e-12,(float)6.63693e-13,(float)
	    -2.15126e-13,(float)7.2659e-14,(float)-2.5465e-14,(float)
	    9.229e-15,(float)-3.448e-15,(float)1.325e-15,(float)-5.22e-16,(
	    float)2.1e-16,(float)-8.7e-17,(float)3.6e-17 };
    static real twodpi = (float).63661977236758134;
    static real pi4 = (float).78539816339744831;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(), log(), sin();

    /* Local variables */
    static real ampl, xmax, xsml;
    extern doublereal besj0_();
    static integer ntth0;
    static real y, z__, theta;
    extern doublereal csevl_();
    extern integer inits_();
    extern doublereal r1mach_();
    extern /* Subroutine */ int xermsg_();
    static integer ntm0, nty0;

/* ***BEGIN PROLOGUE  BESY0 */
/* ***PURPOSE  Compute the Bessel function of the second kind of order */
/*            zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10A1 */
/* ***TYPE      SINGLE PRECISION (BESY0-S, DBESY0-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ZERO, SECOND KIND, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESY0(X) calculates the Bessel function of the second kind */
/* of order zero for real argument X. */

/* Series for BY0        on the interval  0.          to  1.60000D+01 */
/*                                        with weighted error   1.20E-17 
*/
/*                                         log weighted error  16.92 */
/*                               significant figures required  16.15 */
/*                                    decimal places required  17.48 */

/* Series for BM0        on the interval  0.          to  6.25000D-02 */
/*                                        with weighted error   4.98E-17 
*/
/*                                         log weighted error  16.30 */
/*                               significant figures required  14.97 */
/*                                    decimal places required  16.96 */

/* Series for BTH0       on the interval  0.          to  6.25000D-02 */
/*                                        with weighted error   3.67E-17 
*/
/*                                         log weighted error  16.44 */
/*                               significant figures required  15.53 */
/*                                    decimal places required  17.13 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESJ0, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESY0 */
/* ***FIRST EXECUTABLE STATEMENT  BESY0 */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	nty0 = inits_(by0cs, &c__13, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntm0 = inits_(bm0cs, &c__21, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntth0 = inits_(bth0cs, &c__24, &r__1);

	xsml = sqrt(r1mach_(&c__3) * (float)4.);
	xmax = (float)1. / r1mach_(&c__4);
    }
    first = FALSE_;

    if (*x <= (float)0.) {
	xermsg_("SLATEC", "BESY0", "X IS ZERO OR NEGATIVE", &c__1, &c__2, 6L, 
		5L, 21L);
    }
    if (*x > (float)4.) {
	goto L20;
    }

    y = (float)0.;
    if (*x > xsml) {
	y = *x * *x;
    }
    r__1 = y * (float).125 - (float)1.;
    ret_val = twodpi * log(*x * (float).5) * besj0_(x) + (float).375 + csevl_(
	    &r__1, by0cs, &nty0);
    return ret_val;

L20:
    if (*x > xmax) {
	xermsg_("SLATEC", "BESY0", "NO PRECISION BECAUSE X IS BIG", &c__2, &
		c__2, 6L, 5L, 29L);
    }

/* Computing 2nd power */
    r__1 = *x;
    z__ = (float)32. / (r__1 * r__1) - (float)1.;
    ampl = (csevl_(&z__, bm0cs, &ntm0) + (float).75) / sqrt(*x);
    theta = *x - pi4 + csevl_(&z__, bth0cs, &ntth0) / *x;
    ret_val = ampl * sin(theta);

    return ret_val;
} /* besy0_ */

