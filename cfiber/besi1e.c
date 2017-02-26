/* besi1e.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__11 = 11;
static integer c__3 = 3;
static integer c__21 = 21;
static integer c__22 = 22;
static integer c__1 = 1;

/* DECK BESI1E */
doublereal besi1e_(x)
real *x;
{
    /* Initialized data */

    static real bi1cs[11] = { (float)-.001971713261099859,(float)
	    .40734887667546481,(float).034838994299959456,(float)
	    .001545394556300123,(float)4.1888521098377e-5,(float)
	    7.64902676483e-7,(float)1.0042493924e-8,(float)9.9322077e-11,(
	    float)7.6638e-13,(float)4.741e-15,(float)2.4e-17 };
    static real ai1cs[21] = { (float)-.02846744181881479,(float)
	    -.01922953231443221,(float)-6.1151858579437e-4,(float)
	    -2.06997125335e-5,(float)8.58561914581e-6,(float)1.04949824671e-6,
	    (float)-2.9183389184e-7,(float)-1.559378146e-8,(float)
	    1.318012367e-8,(float)-1.44842341e-9,(float)-2.9085122e-10,(float)
	    1.2663889e-10,(float)-1.664947e-11,(float)-1.66665e-12,(float)
	    1.2426e-12,(float)-2.7315e-13,(float)2.023e-14,(float)7.3e-15,(
	    float)-3.33e-15,(float)7.1e-16,(float)-6e-17 };
    static real ai12cs[22] = { (float).02857623501828014,(float)
	    -.00976109749136147,(float)-1.1058893876263e-4,(float)
	    -3.88256480887e-6,(float)-2.5122362377e-7,(float)-2.631468847e-8,(
	    float)-3.83538039e-9,(float)-5.5897433e-10,(float)-1.897495e-11,(
	    float)3.252602e-11,(float)1.41258e-11,(float)2.03564e-12,(float)
	    -7.1985e-13,(float)-4.0836e-13,(float)-2.101e-14,(float)4.273e-14,
	    (float)1.041e-14,(float)-3.82e-15,(float)-1.86e-15,(float)3.3e-16,
	    (float)2.8e-16,(float)-3e-17 };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(), exp(), r_sign();

    /* Local variables */
    static real xmin, xsml;
    static integer ntai1;
    static real y;
    static integer ntai12;
    extern doublereal csevl_();
    extern integer inits_();
    extern doublereal r1mach_();
    extern /* Subroutine */ int xermsg_();
    static integer nti1;

/* ***BEGIN PROLOGUE  BESI1E */
/* ***PURPOSE  Compute the exponentially scaled modified (hyperbolic) */
/*            Bessel function of the first kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESI1E-S, DBSI1E-D) */
/* ***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB, */
/*             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION, */
/*             ORDER ONE, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESI1E(X) calculates the exponentially scaled modified (hyperbolic) */
/* Bessel function of the first kind of order one for real argument X; */
/* i.e., EXP(-ABS(X))*I1(X). */

/* Series for BI1        on the interval  0.          to  9.00000D+00 */
/*                                        with weighted error   2.40E-17 
*/
/*                                         log weighted error  16.62 */
/*                               significant figures required  16.23 */
/*                                    decimal places required  17.14 */

/* Series for AI1        on the interval  1.25000D-01 to  3.33333D-01 */
/*                                        with weighted error   6.98E-17 
*/
/*                                         log weighted error  16.16 */
/*                               significant figures required  14.53 */
/*                                    decimal places required  16.82 */

/* Series for AI12       on the interval  0.          to  1.25000D-01 */
/*                                        with weighted error   3.55E-17 
*/
/*                                         log weighted error  16.45 */
/*                               significant figures required  14.69 */
/*                                    decimal places required  17.12 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890210  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  BESI1E */
/* ***FIRST EXECUTABLE STATEMENT  BESI1E */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	nti1 = inits_(bi1cs, &c__11, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntai1 = inits_(ai1cs, &c__21, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntai12 = inits_(ai12cs, &c__22, &r__1);

	xmin = r1mach_(&c__1) * (float)2.;
	xsml = sqrt(r1mach_(&c__3) * (float)4.5);
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
	xermsg_("SLATEC", "BESI1E", "ABS(X) SO SMALL I1 UNDERFLOWS", &c__1, &
		c__1, 6L, 6L, 29L);
    }
    if (y > xmin) {
	ret_val = *x * (float).5;
    }
    if (y > xsml) {
	r__1 = y * y / (float)4.5 - (float)1.;
	ret_val = *x * (csevl_(&r__1, bi1cs, &nti1) + (float).875);
    }
    ret_val = exp(-y) * ret_val;
    return ret_val;

L20:
    if (y <= (float)8.) {
	r__1 = ((float)48. / y - (float)11.) / (float)5.;
	ret_val = (csevl_(&r__1, ai1cs, &ntai1) + (float).375) / sqrt(y);
    }
    if (y > (float)8.) {
	r__1 = (float)16. / y - (float)1.;
	ret_val = (csevl_(&r__1, ai12cs, &ntai12) + (float).375) / sqrt(y);
    }
    ret_val = r_sign(&ret_val, x);

    return ret_val;
} /* besi1e_ */

