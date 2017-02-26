#include "f2c.h"

/* Table of constant values */

static integer c__23 = 23;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK GAMMA */
doublereal gamma_(x)
real *x;
{
    /* Initialized data */

    static real gcs[23] = { (float).008571195590989331,(float)
	    .004415381324841007,(float).05685043681599363,(float)
	    -.004219835396418561,(float).00132680818121246,(float)
	    -1.89302452979888e-4,(float)3.60692532744124e-5,(float)
	    -6.0567619044608e-6,(float)1.0558295463022e-6,(float)
	    -1.811967365542e-7,(float)3.11772496471e-8,(float)-5.354219639e-9,
	    (float)9.193275519e-10,(float)-1.57794128e-10,(float)
	    2.70798062e-11,(float)-4.6468186e-12,(float)7.97335e-13,(float)
	    -1.368078e-13,(float)2.34731e-14,(float)-4.0274e-15,(float)
	    6.91e-16,(float)-1.185e-16,(float)2.03e-17 };
    static real pi = (float)3.14159265358979324;
    static real sq2pil = (float).91893853320467274;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2;

    /* Builtin functions */
    double sqrt(), r_int(), log(), exp(), sin();

    /* Local variables */
    static integer ngcs;
    static real xmin, xmax;
    static integer i__, n;
    static real y;
    extern doublereal csevl_();
    static real dxrel;
    extern integer inits_();
    extern doublereal r1mach_(), r9lgmc_();
    extern /* Subroutine */ int gamlim_(), xermsg_();
    static real sinpiy;

/* ***BEGIN PROLOGUE  GAMMA */
/* ***PURPOSE  Compute the complete Gamma function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      SINGLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-C) */
/* ***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* GAMMA computes the gamma function at X, where X is not 0, -1, -2, .... 
*/
/* GAMMA and X are single precision. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, GAMLIM, INITS, R1MACH, R9LGMC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  GAMMA */
/* SQ2PIL IS LOG (SQRT (2.*PI) ) */

/* LANL DEPENDENT CODE REMOVED 81.02.04 */

/* ***FIRST EXECUTABLE STATEMENT  GAMMA */
    if (first) {

/* ------------------------------------------------------------------
--- */
/* INITIALIZE.  FIND LEGAL BOUNDS FOR X, AND DETERMINE THE NUMBER OF 
*/
/* TERMS IN THE SERIES REQUIRED TO ATTAIN AN ACCURACY TEN TIMES BETTER
 */
/* THAN MACHINE PRECISION. */

	r__1 = r1mach_(&c__3) * (float).1;
	ngcs = inits_(gcs, &c__23, &r__1);

	gamlim_(&xmin, &xmax);
	dxrel = sqrt(r1mach_(&c__4));

/* ------------------------------------------------------------------
--- */
/* FINISH INITIALIZATION.  START EVALUATING GAMMA(X). */

    }
    first = FALSE_;

    y = dabs(*x);
    if (y > (float)10.) {
	goto L50;
    }

/* COMPUTE GAMMA(X) FOR ABS(X) .LE. 10.0.  REDUCE INTERVAL AND */
/* FIND GAMMA(1+Y) FOR 0. .LE. Y .LT. 1. FIRST OF ALL. */

    n = *x;
    if (*x < (float)0.) {
	--n;
    }
    y = *x - n;
    --n;
    r__1 = y * (float)2. - (float)1.;
    ret_val = csevl_(&r__1, gcs, &ngcs) + (float).9375;
    if (n == 0) {
	return ret_val;
    }

    if (n > 0) {
	goto L30;
    }

/* COMPUTE GAMMA(X) FOR X .LT. 1. */

    n = -n;
    if (*x == (float)0.) {
	xermsg_("SLATEC", "GAMMA", "X IS 0", &c__4, &c__2, 6L, 5L, 6L);
    }
    if (*x < (float)0. && *x + n - 2 == (float)0.) {
	xermsg_("SLATEC", "GAMMA", "X IS A NEGATIVE INTEGER", &c__4, &c__2, 
		6L, 5L, 23L);
    }
    r__2 = *x - (float).5;
    if (*x < (float)-.5 && (r__1 = (*x - r_int(&r__2)) / *x, dabs(r__1)) < 
	    dxrel) {
	xermsg_("SLATEC", "GAMMA", "ANSWER LT HALF PRECISION BECAUSE X TOO N\
EAR NEGATIVE INTEGER", &c__1, &c__1, 6L, 5L, 60L);
    }

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val /= *x + i__ - 1;
/* L20: */
    }
    return ret_val;

/* GAMMA(X) FOR X .GE. 2. */

L30:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val = (y + i__) * ret_val;
/* L40: */
    }
    return ret_val;

/* COMPUTE GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X). */

L50:
    if (*x > xmax) {
	xermsg_("SLATEC", "GAMMA", "X SO BIG GAMMA OVERFLOWS", &c__3, &c__2, 
		6L, 5L, 24L);
    }

    ret_val = (float)0.;
    if (*x < xmin) {
	xermsg_("SLATEC", "GAMMA", "X SO SMALL GAMMA UNDERFLOWS", &c__2, &
		c__1, 6L, 5L, 27L);
    }
    if (*x < xmin) {
	return ret_val;
    }

    ret_val = exp((y - (float).5) * log(y) - y + sq2pil + r9lgmc_(&y));
    if (*x > (float)0.) {
	return ret_val;
    }

    r__2 = *x - (float).5;
    if ((r__1 = (*x - r_int(&r__2)) / *x, dabs(r__1)) < dxrel) {
	xermsg_("SLATEC", "GAMMA", "ANSWER LT HALF PRECISION, X TOO NEAR NEG\
ATIVE INTEGER", &c__1, &c__1, 6L, 5L, 53L);
    }

    sinpiy = sin(pi * y);
    if (sinpiy == (float)0.) {
	xermsg_("SLATEC", "GAMMA", "X IS A NEGATIVE INTEGER", &c__4, &c__2, 
		6L, 5L, 23L);
    }

    ret_val = -pi / (y * sinpiy * ret_val);

    return ret_val;
} /* gamma_ */

