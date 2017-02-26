/* besk1e.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__11 = 11;
static integer c__3 = 3;
static integer c__17 = 17;
static integer c__14 = 14;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK BESK1E */
doublereal besk1e_(x)
real *x;
{
    /* Initialized data */

    static real bk1cs[11] = { (float).0253002273389477705,(float)
	    -.353155960776544876,(float)-.122611180822657148,(float)
	    -.0069757238596398643,(float)-1.730288957513052e-4,(float)
	    -2.4334061415659e-6,(float)-2.21338763073e-8,(float)
	    -1.411488392e-10,(float)-6.666901e-13,(float)-2.4274e-15,(float)
	    -7e-18 };
    static real ak1cs[17] = { (float).2744313406973883,(float)
	    .07571989953199368,(float)-.0014410515564754,(float)
	    6.650116955125e-5,(float)-4.36998470952e-6,(float)3.5402774997e-7,
	    (float)-3.311163779e-8,(float)3.44597758e-9,(float)-3.8989323e-10,
	    (float)4.720819e-11,(float)-6.04783e-12,(float)8.1284e-13,(float)
	    -1.1386e-13,(float)1.654e-14,(float)-2.48e-15,(float)3.8e-16,(
	    float)-6e-17 };
    static real ak12cs[14] = { (float).06379308343739001,(float)
	    .02832887813049721,(float)-2.4753706739052e-4,(float)
	    5.7719724516e-6,(float)-2.0689392195e-7,(float)9.73998344e-9,(
	    float)-5.5853361e-10,(float)3.732996e-11,(float)-2.82505e-12,(
	    float)2.372e-13,(float)-2.176e-14,(float)2.15e-15,(float)-2.2e-16,
	    (float)2e-17 };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double log(), exp(), sqrt();

    /* Local variables */
    static real xmin, xsml;
    extern doublereal besi1_();
    static integer ntak1;
    static real y;
    static integer ntak12;
    extern doublereal csevl_();
    extern integer inits_();
    extern doublereal r1mach_();
    extern /* Subroutine */ int xermsg_();
    static integer ntk1;

/* ***BEGIN PROLOGUE  BESK1E */
/* ***PURPOSE  Compute the exponentially scaled modified (hyperbolic) */
/*            Bessel function of the third kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESK1E-S, DBSK1E-D) */
/* ***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION, 
*/
/*             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESK1E(X) computes the exponentially scaled modified (hyperbolic) */
/* Bessel function of third kind of order one for real argument */
/* X .GT. 0.0, i.e., EXP(X)*K1(X). */

/* Series for BK1        on the interval  0.          to  4.00000D+00 */
/*                                        with weighted error   7.02E-18 
*/
/*                                         log weighted error  17.15 */
/*                               significant figures required  16.73 */
/*                                    decimal places required  17.67 */

/* Series for AK1        on the interval  1.25000D-01 to  5.00000D-01 */
/*                                        with weighted error   6.06E-17 
*/
/*                                         log weighted error  16.22 */
/*                               significant figures required  15.41 */
/*                                    decimal places required  16.83 */

/* Series for AK12       on the interval  0.          to  1.25000D-01 */
/*                                        with weighted error   2.58E-17 
*/
/*                                         log weighted error  16.59 */
/*                               significant figures required  15.22 */
/*                                    decimal places required  17.16 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI1, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESK1E */
/* ***FIRST EXECUTABLE STATEMENT  BESK1E */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	ntk1 = inits_(bk1cs, &c__11, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntak1 = inits_(ak1cs, &c__17, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntak12 = inits_(ak12cs, &c__14, &r__1);

/* Computing MAX */
	r__1 = log(r1mach_(&c__1)), r__2 = -log(r1mach_(&c__2));
	xmin = exp(dmax(r__1,r__2) + (float).01);
	xsml = sqrt(r1mach_(&c__3) * (float)4.);
    }
    first = FALSE_;

    if (*x <= (float)0.) {
	xermsg_("SLATEC", "BESK1E", "X IS ZERO OR NEGATIVE", &c__2, &c__2, 6L,
		 6L, 21L);
    }
    if (*x > (float)2.) {
	goto L20;
    }

    if (*x < xmin) {
	xermsg_("SLATEC", "BESK1E", "X SO SMALL K1 OVERFLOWS", &c__3, &c__2, 
		6L, 6L, 23L);
    }
    y = (float)0.;
    if (*x > xsml) {
	y = *x * *x;
    }
    r__1 = y * (float).5 - (float)1.;
    ret_val = exp(*x) * (log(*x * (float).5) * besi1_(x) + (csevl_(&r__1, 
	    bk1cs, &ntk1) + (float).75) / *x);
    return ret_val;

L20:
    if (*x <= (float)8.) {
	r__1 = ((float)16. / *x - (float)5.) / (float)3.;
	ret_val = (csevl_(&r__1, ak1cs, &ntak1) + (float)1.25) / sqrt(*x);
    }
    if (*x > (float)8.) {
	r__1 = (float)16. / *x - (float)1.;
	ret_val = (csevl_(&r__1, ak12cs, &ntak12) + (float)1.25) / sqrt(*x);
    }

    return ret_val;
} /* besk1e_ */

