/* besi0e.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__12 = 12;
static integer c__3 = 3;
static integer c__21 = 21;
static integer c__22 = 22;

/* DECK BESI0E */
doublereal besi0e_(x)
real *x;
{
    /* Initialized data */

    static real bi0cs[12] = { (float)-.07660547252839144951,(float)
	    1.92733795399380827,(float).2282644586920301339,(float)
	    .01304891466707290428,(float)4.3442709008164874e-4,(float)
	    9.42265768600193e-6,(float)1.4340062895106e-7,(float)
	    1.61384906966e-9,(float)1.396650044e-11,(float)9.579451e-14,(
	    float)5.3339e-16,(float)2.45e-18 };
    static real ai0cs[21] = { (float).07575994494023796,(float)
	    .00759138081082334,(float)4.1531313389237e-4,(float)
	    1.070076463439e-5,(float)-7.90117997921e-6,(float)
	    -7.8261435014e-7,(float)2.7838499429e-7,(float)8.2524726e-9,(
	    float)-1.204463945e-8,(float)1.55964859e-9,(float)2.2925563e-10,(
	    float)-1.1916228e-10,(float)1.757854e-11,(float)1.12822e-12,(
	    float)-1.14684e-12,(float)2.7155e-13,(float)-2.415e-14,(float)
	    -6.08e-15,(float)3.14e-15,(float)-7.1e-16,(float)7e-17 };
    static real ai02cs[22] = { (float).05449041101410882,(float)
	    .00336911647825569,(float)6.889758346918e-5,(float)
	    2.89137052082e-6,(float)2.0489185893e-7,(float)2.266668991e-8,(
	    float)3.39623203e-9,(float)4.9406022e-10,(float)1.188914e-11,(
	    float)-3.149915e-11,(float)-1.32158e-11,(float)-1.79419e-12,(
	    float)7.1801e-13,(float)3.8529e-13,(float)1.539e-14,(float)
	    -4.151e-14,(float)-9.54e-15,(float)3.82e-15,(float)1.76e-15,(
	    float)-3.4e-16,(float)-2.7e-16,(float)3e-17 };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(), exp();

    /* Local variables */
    static real xsml;
    static integer ntai0;
    static real y;
    static integer ntai02;
    extern doublereal csevl_();
    extern integer inits_();
    extern doublereal r1mach_();
    static integer nti0;

/* ***BEGIN PROLOGUE  BESI0E */
/* ***PURPOSE  Compute the exponentially scaled modified (hyperbolic) */
/*            Bessel function of the first kind of order zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESI0E-S, DBSI0E-D) */
/* ***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB, */
/*             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION, */
/*             ORDER ZERO, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESI0E(X) calculates the exponentially scaled modified (hyperbolic) */
/* Bessel function of the first kind of order zero for real argument X; */
/* i.e., EXP(-ABS(X))*I0(X). */


/* Series for BI0        on the interval  0.          to  9.00000D+00 */
/*                                        with weighted error   2.46E-18 
*/
/*                                         log weighted error  17.61 */
/*                               significant figures required  17.90 */
/*                                    decimal places required  18.15 */


/* Series for AI0        on the interval  1.25000D-01 to  3.33333D-01 */
/*                                        with weighted error   7.87E-17 
*/
/*                                         log weighted error  16.10 */
/*                               significant figures required  14.69 */
/*                                    decimal places required  16.76 */


/* Series for AI02       on the interval  0.          to  1.25000D-01 */
/*                                        with weighted error   3.79E-17 
*/
/*                                         log weighted error  16.42 */
/*                               significant figures required  14.86 */
/*                                    decimal places required  17.09 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890313  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  BESI0E */
/* ***FIRST EXECUTABLE STATEMENT  BESI0E */
    if (first) {
	r__1 = r1mach_(&c__3) * (float).1;
	nti0 = inits_(bi0cs, &c__12, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntai0 = inits_(ai0cs, &c__21, &r__1);
	r__1 = r1mach_(&c__3) * (float).1;
	ntai02 = inits_(ai02cs, &c__22, &r__1);
	xsml = sqrt(r1mach_(&c__3) * (float)4.5);
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > (float)3.) {
	goto L20;
    }

    ret_val = (float)1. - *x;
    if (y > xsml) {
	r__1 = y * y / (float)4.5 - (float)1.;
	ret_val = exp(-y) * (csevl_(&r__1, bi0cs, &nti0) + (float)2.75);
    }
    return ret_val;

L20:
    if (y <= (float)8.) {
	r__1 = ((float)48. / y - (float)11.) / (float)5.;
	ret_val = (csevl_(&r__1, ai0cs, &ntai0) + (float).375) / sqrt(y);
    }
    if (y > (float)8.) {
	r__1 = (float)16. / y - (float)1.;
	ret_val = (csevl_(&r__1, ai02cs, &ntai02) + (float).375) / sqrt(y);
    }

    return ret_val;
} /* besi0e_ */

