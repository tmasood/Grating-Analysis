/* besj0.f -- translated by f2c (version 19961209).
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

/* DECK BESJ0 */
doublereal besj0_(x)
real *x;
{
    /* Initialized data */

    static real bj0cs[13] = { (float).100254161968939137,(float)
	    -.665223007764405132,(float).248983703498281314,(float)
	    -.0332527231700357697,(float).0023114179304694015,(float)
	    -9.9112774199508e-5,(float)2.8916708643998e-6,(float)
	    -6.1210858663e-8,(float)9.838650793e-10,(float)-1.24235515e-11,(
	    float)1.265433e-13,(float)-1.0619e-15,(float)7.4e-18 };
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
    static real pi4 = (float).78539816339744831;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(), cos();

    /* Local variables */
    static real ampl, xmax, xsml;
    static integer ntth0;
    static real y, z__, theta;
    extern doublereal csevl_();
    extern integer inits_();
    extern doublereal r1mach_();
    extern /* Subroutine */ int xermsg_();
    static integer ntj0, ntm0;

/* ***BEGIN PROLOGUE  BESJ0 */
/* ***PURPOSE  Compute the Bessel function of the first kind of order */
/*            zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10A1 */
/* ***TYPE      SINGLE PRECISION (BESJ0-S, DBESJ0-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ZERO, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESJ0(X) calculates the Bessel function of the first kind of */
/* order zero for real argument X. */

/* Series for BJ0        on the interval  0.          to  1.60000D+01 */
/*                                        with weighted error   7.47E-18 
*/
/*                                         log weighted error  17.13 */
/*                               significant figures required  16.98 */
/*                                    decimal places required  17.68 */

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
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890210  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESJ0 */
/* ***FIRST EXECUTABLE STATEMENT  BESJ0 */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	ntj0 = inits_(bj0cs, &c__13, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntm0 = inits_(bm0cs, &c__21, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntth0 = inits_(bth0cs, &c__24, &r__1);

	xsml = sqrt(r1mach_(&c__3) * (float)8.);
	xmax = (float)1. / r1mach_(&c__4);
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > (float)4.) {
	goto L20;
    }

    ret_val = (float)1.;
    if (y > xsml) {
	r__1 = y * (float).125 * y - (float)1.;
	ret_val = csevl_(&r__1, bj0cs, &ntj0);
    }
    return ret_val;

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "BESJ0", "NO PRECISION BECAUSE ABS(X) IS TOO BIG", &
		c__1, &c__2, 6L, 5L, 38L);
    }

/* Computing 2nd power */
    r__1 = y;
    z__ = (float)32. / (r__1 * r__1) - (float)1.;
    ampl = (csevl_(&z__, bm0cs, &ntm0) + (float).75) / sqrt(y);
    theta = y - pi4 + csevl_(&z__, bth0cs, &ntth0) / y;
    ret_val = ampl * cos(theta);

    return ret_val;
} /* besj0_ */

