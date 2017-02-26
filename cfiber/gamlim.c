#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* DECK GAMLIM */
/* Subroutine */ int gamlim_(xmin, xmax)
real *xmin, *xmax;
{
    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    double log();

    /* Local variables */
    static real xold;
    static integer i__;
    extern doublereal r1mach_();
    static real alnbig, alnsml;
    extern /* Subroutine */ int xermsg_();
    static real xln;

/* ***BEGIN PROLOGUE  GAMLIM */
/* ***PURPOSE  Compute the minimum and maximum bounds for the argument in 
*/
/*            the Gamma function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A, R2 */
/* ***TYPE      SINGLE PRECISION (GAMLIM-S, DGAMLM-D) */
/* ***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS 
*/
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Calculate the minimum and maximum legal bounds for X in GAMMA(X). */
/* XMIN and XMAX are not the only bounds, but they are the only non- */
/* trivial ones to calculate. */

/*             Output Arguments -- */
/* XMIN   minimum legal value of X in GAMMA(X).  Any smaller value of */
/*        X might result in underflow. */
/* XMAX   maximum legal value of X in GAMMA(X).  Any larger value will */
/*        cause overflow. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  GAMLIM */
/* ***FIRST EXECUTABLE STATEMENT  GAMLIM */
    alnsml = log(r1mach_(&c__1));
    *xmin = -alnsml;
    for (i__ = 1; i__ <= 10; ++i__) {
	xold = *xmin;
	xln = log(*xmin);
	*xmin -= *xmin * ((*xmin + (float).5) * xln - *xmin - (float).2258 + 
		alnsml) / (*xmin * xln + (float).5);
	if ((r__1 = *xmin - xold, dabs(r__1)) < (float).005) {
	    goto L20;
	}
/* L10: */
    }
    xermsg_("SLATEC", "GAMLIM", "UNABLE TO FIND XMIN", &c__1, &c__2, 6L, 6L, 
	    19L);

L20:
    *xmin = -(*xmin) + (float).01;

    alnbig = log(r1mach_(&c__2));
    *xmax = alnbig;
    for (i__ = 1; i__ <= 10; ++i__) {
	xold = *xmax;
	xln = log(*xmax);
	*xmax -= *xmax * ((*xmax - (float).5) * xln - *xmax + (float).9189 - 
		alnbig) / (*xmax * xln - (float).5);
	if ((r__1 = *xmax - xold, dabs(r__1)) < (float).005) {
	    goto L40;
	}
/* L30: */
    }
    xermsg_("SLATEC", "GAMLIM", "UNABLE TO FIND XMAX", &c__2, &c__2, 6L, 6L, 
	    19L);

L40:
    *xmax += (float)-.01;
/* Computing MAX */
    r__1 = *xmin, r__2 = -(*xmax) + (float)1.;
    *xmin = dmax(r__1,r__2);

    return 0;
} /* gamlim_ */

