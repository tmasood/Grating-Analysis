/* besj1.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__12 = 12;
static integer c__3 = 3;
static integer c__21 = 21;
static integer c__24 = 24;
static integer c__1 = 1;
static integer c__4 = 4;
static integer c__2 = 2;

/* DECK BESJ1 */
doublereal besj1_(x)
real *x;
{
    /* Initialized data */

    static real bj1cs[12] = { (float)-.11726141513332787,(float)
	    -.2536152183079064,(float).050127080984469569,(float)
	    -.004631514809625081,(float)2.47996229415914e-4,(float)
	    -8.678948686278e-6,(float)2.14293917143e-7,(float)-3.936093079e-9,
	    (float)5.5911823e-11,(float)-6.32761e-13,(float)5.84e-15,(float)
	    -4.4e-17 };
    static real bm1cs[21] = { (float).1047362510931285,(float)
	    .00442443893702345,(float)-5.661639504035e-5,(float)
	    2.31349417339e-6,(float)-1.7377182007e-7,(float)1.89320993e-8,(
	    float)-2.65416023e-9,(float)4.4740209e-10,(float)-8.691795e-11,(
	    float)1.891492e-11,(float)-4.51884e-12,(float)1.16765e-12,(float)
	    -3.2265e-13,(float)9.45e-14,(float)-2.913e-14,(float)9.39e-15,(
	    float)-3.15e-15,(float)1.09e-15,(float)-3.9e-16,(float)1.4e-16,(
	    float)-5e-17 };
    static real bth1cs[24] = { (float).7406014102631385,(float)
	    -.00457175565963769,(float)1.19818510964326e-4,(float)
	    -6.964561891648e-6,(float)6.55495621447e-7,(float)
	    -8.4066228945e-8,(float)1.3376886564e-8,(float)-2.499565654e-9,(
	    float)5.294951e-10,(float)-1.24135944e-10,(float)3.1656485e-11,(
	    float)-8.66864e-12,(float)2.523758e-12,(float)-7.75085e-13,(float)
	    2.49527e-13,(float)-8.3773e-14,(float)2.9205e-14,(float)
	    -1.0534e-14,(float)3.919e-15,(float)-1.5e-15,(float)5.89e-16,(
	    float)-2.37e-16,(float)9.7e-17,(float)-4e-17 };
    static real pi4 = (float).78539816339744831;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(), r_sign(), cos();

    /* Local variables */
    static real ampl, xmin, xmax, xsml;
    static integer ntth1;
    static real y, z__, theta;
    extern doublereal csevl_();
    extern integer inits_();
    extern doublereal r1mach_();
    extern /* Subroutine */ int xermsg_();
    static integer ntj1, ntm1;

/* ***BEGIN PROLOGUE  BESJ1 */
/* ***PURPOSE  Compute the Bessel function of the first kind of order one.
 */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10A1 */
/* ***TYPE      SINGLE PRECISION (BESJ1-S, DBESJ1-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ONE, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESJ1(X) calculates the Bessel function of the first kind of */
/* order one for real argument X. */

/* Series for BJ1        on the interval  0.          to  1.60000D+01 */
/*                                        with weighted error   4.48E-17 
*/
/*                                         log weighted error  16.35 */
/*                               significant figures required  15.77 */
/*                                    decimal places required  16.89 */

/* Series for BM1        on the interval  0.          to  6.25000D-02 */
/*                                        with weighted error   5.61E-17 
*/
/*                                         log weighted error  16.25 */
/*                               significant figures required  14.97 */
/*                                    decimal places required  16.91 */

/* Series for BTH1       on the interval  0.          to  6.25000D-02 */
/*                                        with weighted error   4.10E-17 
*/
/*                                         log weighted error  16.39 */
/*                               significant figures required  15.96 */
/*                                    decimal places required  17.08 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780601  DATE WRITTEN */
/*   890210  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESJ1 */
/* ***FIRST EXECUTABLE STATEMENT  BESJ1 */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	ntj1 = inits_(bj1cs, &c__12, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntm1 = inits_(bm1cs, &c__21, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntth1 = inits_(bth1cs, &c__24, &r__1);

	xsml = sqrt(r1mach_(&c__3) * (float)8.);
	xmin = r1mach_(&c__1) * (float)2.;
	xmax = (float)1. / r1mach_(&c__4);
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > (float)4.) {
	goto L20;
    }

    ret_val = (float)0.;
    if (y == (float)0.) {
	return ret_val;
    }
    if (y <= xmin) {
	xermsg_("SLATEC", "BESJ1", "ABS(X) SO SMALL J1 UNDERFLOWS", &c__1, &
		c__1, 6L, 5L, 29L);
    }
    if (y > xmin) {
	ret_val = *x * (float).5;
    }
    if (y > xsml) {
	r__1 = y * (float).125 * y - (float)1.;
	ret_val = *x * (csevl_(&r__1, bj1cs, &ntj1) + (float).25);
    }
    return ret_val;

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "BESJ1", "NO PRECISION BECAUSE ABS(X) IS TOO BIG", &
		c__2, &c__2, 6L, 5L, 38L);
    }
/* Computing 2nd power */
    r__1 = y;
    z__ = (float)32. / (r__1 * r__1) - (float)1.;
    ampl = (csevl_(&z__, bm1cs, &ntm1) + (float).75) / sqrt(y);
    theta = y - pi4 * (float)3. + csevl_(&z__, bth1cs, &ntth1) / y;
    ret_val = r_sign(&ampl, x) * cos(theta);

    return ret_val;
} /* besj1_ */

