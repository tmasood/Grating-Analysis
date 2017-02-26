/* besy1.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__14 = 14;
static integer c__3 = 3;
static integer c__21 = 21;
static integer c__24 = 24;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__4 = 4;

/* DECK BESY1 */
doublereal besy1_(x)
real *x;
{
    /* Initialized data */

    static real by1cs[14] = { (float).03208047100611908629,(float)
	    1.26270789743350045,(float).006499961899923175,(float)
	    -.08936164528860504117,(float).01325088122175709545,(float)
	    -8.9790591196483523e-4,(float)3.647361487958306e-5,(float)
	    -1.001374381666e-6,(float)1.99453965739e-8,(float)
	    -3.0230656018e-10,(float)3.60987815e-12,(float)-3.487488e-14,(
	    float)2.7838e-16,(float)-1.86e-18 };
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
    static real twodpi = (float).63661977236758134;
    static real pi4 = (float).78539816339744831;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double log(), exp(), sqrt(), sin();

    /* Local variables */
    static real ampl, xmin, xmax, xsml;
    extern doublereal besj1_();
    static integer ntth1;
    static real y, z__, theta;
    extern doublereal csevl_();
    extern integer inits_();
    extern doublereal r1mach_();
    extern /* Subroutine */ int xermsg_();
    static integer ntm1, nty1;

/* ***BEGIN PROLOGUE  BESY1 */
/* ***PURPOSE  Compute the Bessel function of the second kind of order */
/*            one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10A1 */
/* ***TYPE      SINGLE PRECISION (BESY1-S, DBESY1-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ONE, SECOND KIND, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESY1(X) calculates the Bessel function of the second kind of */
/* order one for real argument X. */

/* Series for BY1        on the interval  0.          to  1.60000D+01 */
/*                                        with weighted error   1.87E-18 
*/
/*                                         log weighted error  17.73 */
/*                               significant figures required  17.83 */
/*                                    decimal places required  18.30 */

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
/* ***ROUTINES CALLED  BESJ1, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESY1 */
/* ***FIRST EXECUTABLE STATEMENT  BESY1 */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	nty1 = inits_(by1cs, &c__14, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntm1 = inits_(bm1cs, &c__21, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntth1 = inits_(bth1cs, &c__24, &r__1);

/* Computing MAX */
	r__1 = log(r1mach_(&c__1)), r__2 = -log(r1mach_(&c__2));
	xmin = exp(dmax(r__1,r__2) + (float).01) * (float)1.571;
	xsml = sqrt(r1mach_(&c__3) * (float)4.);
	xmax = (float)1. / r1mach_(&c__4);
    }
    first = FALSE_;

    if (*x <= (float)0.) {
	xermsg_("SLATEC", "BESY1", "X IS ZERO OR NEGATIVE", &c__1, &c__2, 6L, 
		5L, 21L);
    }
    if (*x > (float)4.) {
	goto L20;
    }

    if (*x < xmin) {
	xermsg_("SLATEC", "BESY1", "X SO SMALL Y1 OVERFLOWS", &c__3, &c__2, 
		6L, 5L, 23L);
    }
    y = (float)0.;
    if (*x > xsml) {
	y = *x * *x;
    }
    r__1 = y * (float).125 - (float)1.;
    ret_val = twodpi * log(*x * (float).5) * besj1_(x) + (csevl_(&r__1, by1cs,
	     &nty1) + (float).5) / *x;
    return ret_val;

L20:
    if (*x > xmax) {
	xermsg_("SLATEC", "BESY1", "NO PRECISION BECAUSE X IS BIG", &c__2, &
		c__2, 6L, 5L, 29L);
    }

/* Computing 2nd power */
    r__1 = *x;
    z__ = (float)32. / (r__1 * r__1) - (float)1.;
    ampl = (csevl_(&z__, bm1cs, &ntm1) + (float).75) / sqrt(*x);
    theta = *x - pi4 * (float)3. + csevl_(&z__, bth1cs, &ntth1) / *x;
    ret_val = ampl * sin(theta);

    return ret_val;
} /* besy1_ */

