#include<stdio.h>
#include<math.h>
#include "f2c.h"

/* Table of constant values */

static integer c__15 = 15;
static integer c__5 = 5;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__6 = 6;

/* DECK BESY */
/* Subroutine */ int besy(real *x, real *fnu, int *n, real *y)
{
    /* Initialized data */

    static integer nulim[2] = { 70,100 };

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(), log();

    /* Local variables */
    static real elim;
    static integer iflw;
    static real xlim;
    extern doublereal besy0_(), besy1_();
    static integer i__, j;
    static real s, w[2], flgjy;
    extern /* Subroutine */ int yairy_();
    extern /* Subroutine */ int asyjy_();
    extern integer i1mach_();
    static real s1, s2;
    extern doublereal r1mach_();
    static integer nb;
    static real cn;
    static integer nd;
    static real fn;
    static integer nn;
    static real tm, wk[7];
    extern /* Subroutine */ int besynu_(), xermsg_();
    static real w2n, ran;
    static integer nud;
    static real dnu, azn, trx, xxn;

/* ***BEGIN PROLOGUE  BESY */
/* ***PURPOSE  Implement forward recursion on the three term recursion */
/*            relation for a sequence of non-negative order Bessel */
/*            functions Y/SUB(FNU+I-1)/(X), I=1,...,N for real, positive 
*/
/*            X and non-negative orders FNU. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10A3 */
/* ***TYPE      SINGLE PRECISION (BESY-S, DBESY-D) */
/* ***KEYWORDS  SPECIAL FUNCTIONS, Y BESSEL FUNCTION */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         BESY implements forward recursion on the three term */
/*         recursion relation for a sequence of non-negative order Bessel 
*/
/*         functions Y/sub(FNU+I-1)/(X), I=1,N for real X .GT. 0.0E0 and 
*/
/*         non-negative orders FNU.  If FNU .LT. NULIM, orders FNU and */
/*         FNU+1 are obtained from BESYNU which computes by a power */
/*         series for X .LE. 2, the K Bessel function of an imaginary */
/*         argument for 2 .LT. X .LE. 20 and the asymptotic expansion for 
*/
/*         X .GT. 20. */

/*         If FNU .GE. NULIM, the uniform asymptotic expansion is coded */
/*         in ASYJY for orders FNU and FNU+1 to start the recursion. */
/*         NULIM is 70 or 100 depending on whether N=1 or N .GE. 2.  An */
/*         overflow test is made on the leading term of the asymptotic */
/*         expansion before any extensive computation is done. */

/*     Description of Arguments */

/*         Input */
/*           X      - X .GT. 0.0E0 */
/*           FNU    - order of the initial Y function, FNU .GE. 0.0E0 */
/*           N      - number of members in the sequence, N .GE. 1 */

/*         Output */
/*           Y      - a vector whose first N components contain values */
/*                    for the sequence Y(I)=Y/sub(FNU+I-1)/(X), I=1,N. */

/*     Error Conditions */
/*         Improper input arguments - a fatal error */
/*         Overflow - a fatal error */

/* ***REFERENCES  F. W. J. Olver, Tables of Bessel Functions of Moderate 
*/
/*                 or Large Orders, NPL Mathematical Tables 6, Her */
/*                 Majesty's Stationery Office, London, 1962. */
/*               N. M. Temme, On the numerical evaluation of the modified 
*/
/*                 Bessel function of the third kind, Journal of */
/*                 Computational Physics 19, (1975), pp. 324-337. */
/*               N. M. Temme, On the numerical evaluation of the ordinary 
*/
/*                 Bessel function of the second kind, Journal of */
/*                 Computational Physics 21, (1976), pp. 343-350. */
/* ***ROUTINES CALLED  ASYJY, BESY0, BESY1, BESYNU, I1MACH, R1MACH, */
/*                    XERMSG, YAIRY */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BESY */

    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  BESY */
    nn = -i1mach_(&c__15);
    elim = (nn * r1mach_(&c__5) - (float)3.) * (float)2.303;
    xlim = r1mach_(&c__1) * (float)1e3;
    if (*fnu < (float)0.) {
	goto L140;
    }
    if (*x <= (float)0.) {
	goto L150;
    }
    if (*x < xlim) {
	goto L170;
    }
    if (*n < 1) {
	goto L160;
    }

/*     ND IS A DUMMY VARIABLE FOR N */

    nd = *n;
    nud = (integer) (*fnu);
    dnu = *fnu - nud;
    nn = min(2,nd);
    fn = *fnu + *n - 1;
    if (fn < (float)2.) {
	goto L100;
    }

/*     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION) */
/*     FOR THE LAST ORDER, FNU+N-1.GE.NULIM */

    xxn = *x / fn;
    w2n = (float)1. - xxn * xxn;
    if (w2n <= (float)0.) {
	goto L10;
    }
    ran = sqrt(w2n);
    azn = log((ran + (float)1.) / xxn) - ran;
    cn = fn * azn;
    if (cn > elim) {
      printf("error elim %f fnu %f n %d x %f\n",elim,*fnu,*n, *x);
      goto L170;
    }
L10:
    if (nud < nulim[nn - 1]) {
	goto L20;
    }

/*     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1.GE.NULIM */

    flgjy = (float)-1.;
    asyjy_(yairy_, x, fnu, &flgjy, &nn, &y[1], wk, &iflw);
    if (iflw != 0) {
	goto L170;
    }
    if (nn == 1) {
	return 0;
    }
    trx = (float)2. / *x;
    tm = (*fnu + *fnu + (float)2.) / *x;
    goto L80;

L20:
    if (dnu != (float)0.) {
	goto L30;
    }
    s1 = besy0_(x);
    if (nud == 0 && nd == 1) {
	goto L70;
    }
    s2 = besy1_(x);
    goto L40;
L30:
    nb = 2;
    if (nud == 0 && nd == 1) {
	nb = 1;
    }
    besynu_(x, &dnu, &nb, w);
    s1 = w[0];
    if (nb == 1) {
	goto L70;
    }
    s2 = w[1];
L40:
    trx = (float)2. / *x;
    tm = (dnu + dnu + (float)2.) / *x;
/*     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2) */
    if (nd == 1) {
	--nud;
    }
    if (nud > 0) {
	goto L50;
    }
    if (nd > 1) {
	goto L70;
    }
    s1 = s2;
    goto L70;
L50:
    i__1 = nud;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = s2;
	s2 = tm * s2 - s1;
	s1 = s;
	tm += trx;
/* L60: */
    }
    if (nd == 1) {
	s1 = s2;
    }
L70:
    y[1] = s1;
    if (nd == 1) {
	return 0;
    }
    y[2] = s2;
L80:
    if (nd == 2) {
	return 0;
    }
/*     FORWARD RECUR FROM FNU+2 TO FNU+N-1 */
    i__1 = nd;
    for (i__ = 3; i__ <= i__1; ++i__) {
	y[i__] = tm * y[i__ - 1] - y[i__ - 2];
	tm += trx;
/* L90: */
    }
    return 0;

L100:
/*     OVERFLOW TEST */
    if (fn <= (float)1.) {
	goto L110;
    }
    if (-fn * (log(*x) - (float).693) > elim) {
	goto L170;
    }
L110:
    if (dnu == (float)0.) {
	goto L120;
    }
    besynu_(x, fnu, &nd, &y[1]);
    return 0;
L120:
    j = nud;
    if (j == 1) {
	goto L130;
    }
    ++j;
    y[j] = besy0_(x);
    if (nd == 1) {
	return 0;
    }
    ++j;
L130:
    y[j] = besy1_(x);
    if (nd == 1) {
	return 0;
    }
    trx = (float)2. / *x;
    tm = trx;
    goto L80;



L140:
    xermsg_("SLATEC", "BESY", "ORDER, FNU, LESS THAN ZERO", &c__2, &c__1, 6L, 
	    4L, 26L);
    return 0;
L150:
    xermsg_("SLATEC", "BESY", "X LESS THAN OR EQUAL TO ZERO", &c__2, &c__1, 
	    6L, 4L, 28L);
    return 0;
L160:
    xermsg_("SLATEC", "BESY", "N LESS THAN ONE", &c__2, &c__1, 6L, 4L, 15L);
    return 0;
L170:
    xermsg_("SLATEC", "BESY", "OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL", &
	    c__6, &c__1, 6L, 4L, 43L);
    return 0;
} /* besy_ */

