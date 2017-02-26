/* r9lgmc.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__6 = 6;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK R9LGMC */
doublereal r9lgmc_(x)
real *x;
{
    /* Initialized data */

    static real algmcs[6] = { (float).166638948045186,(float)
	    -1.38494817606e-5,(float)9.8108256e-9,(float)-1.80912e-11,(float)
	    6.22e-14,(float)-3e-16 };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double sqrt(), log(), exp();

    /* Local variables */
    static real xbig, xmax;
    static integer nalgm;
    extern doublereal csevl_();
    extern integer inits_();
    extern doublereal r1mach_();
    extern /* Subroutine */ int xermsg_();

/* ***BEGIN PROLOGUE  R9LGMC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the log Gamma correction factor so that */
/*            LOG(GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X */
/*            + R9LGMC(X). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7E */
/* ***TYPE      SINGLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C) */
/* ***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB, */
/*             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Compute the log gamma correction factor for X .GE. 10.0 so that */
/*  LOG (GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + R9LGMC(X) */

/* Series for ALGM       on the interval  0.          to  1.00000D-02 */
/*                                        with weighted error   3.40E-16 
*/
/*                                         log weighted error  15.47 */
/*                               significant figures required  14.39 */
/*                                    decimal places required  15.86 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  R9LGMC */
/* ***FIRST EXECUTABLE STATEMENT  R9LGMC */
    if (first) {
	r__1 = r1mach_(&c__3);
	nalgm = inits_(algmcs, &c__6, &r__1);
	xbig = (float)1. / sqrt(r1mach_(&c__3));
/* Computing MIN */
	r__1 = log(r1mach_(&c__2) / (float)12.), r__2 = -log(r1mach_(&c__1) * 
		(float)12.);
	xmax = exp((dmin(r__1,r__2)));
    }
    first = FALSE_;

    if (*x < (float)10.) {
	xermsg_("SLATEC", "R9LGMC", "X MUST BE GE 10", &c__1, &c__2, 6L, 6L, 
		15L);
    }
    if (*x >= xmax) {
	goto L20;
    }

    ret_val = (float)1. / (*x * (float)12.);
    if (*x < xbig) {
/* Computing 2nd power */
	r__2 = (float)10. / *x;
	r__1 = r__2 * r__2 * (float)2. - (float)1.;
	ret_val = csevl_(&r__1, algmcs, &nalgm) / *x;
    }
    return ret_val;

L20:
    ret_val = (float)0.;
    xermsg_("SLATEC", "R9LGMC", "X SO BIG R9LGMC UNDERFLOWS", &c__2, &c__1, 
	    6L, 6L, 26L);
    return ret_val;

} /* r9lgmc_ */

