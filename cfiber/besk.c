#include<stdio.h>
/* besk.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__15 = 15;
/* static integer c__12 = 12; */
static integer c__5 = 5;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__6 = 6;

/* DECK BESK */
/* Subroutine */ int besk(real *x, real *fnu, int *kode, int *n, real *y, int *nz)
{
    /* Initialized data */

    static integer nulim[2] = { 35,70 };

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(), log();

    /* Local variables */
    static real elim, xlim;
    extern doublereal besk0_(), besk1_();
    static integer i__, j, k;
    static real s, t, w[2], flgik;
    extern /* Subroutine */ int asyik_();
    extern integer i1mach_();
    static real s1, s2;
    extern doublereal besk0e_(), besk1e_(), r1mach_();
    static integer nb;
    static real cn;
    static integer nd;
    static real fn;
    static integer nn;
    static real tm;
    static integer mz;
    static real zn;
    extern /* Subroutine */ int besknu_(), xermsg_();
    static real gln, fnn;
    static integer nud;
    static real dnu, gnu, etx, trx, rtz;

/* ***BEGIN PROLOGUE  BESK */
/* ***PURPOSE  Implement forward recursion on the three term recursion */
/*            relation for a sequence of non-negative order Bessel */
/*            functions K/SUB(FNU+I-1)/(X), or scaled Bessel functions */
/*            EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N for real, positive */
/*            X and non-negative orders FNU. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10B3 */
/* ***TYPE      SINGLE PRECISION (BESK-S, DBESK-D) */
/* ***KEYWORDS  K BESSEL FUNCTION, SPECIAL FUNCTIONS */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         BESK implements forward recursion on the three term */
/*         recursion relation for a sequence of non-negative order Bessel 
*/
/*         functions K/sub(FNU+I-1)/(X), or scaled Bessel functions */
/*         EXP(X)*K/sub(FNU+I-1)/(X), I=1,...,N for real X .GT. 0.0E0 and 
*/
/*         non-negative orders FNU.  If FNU .LT. NULIM, orders FNU and */
/*         FNU+1 are obtained from BESKNU to start the recursion.  If */
/*         FNU .GE. NULIM, the uniform asymptotic expansion is used for */
/*         orders FNU and FNU+1 to start the recursion.  NULIM is 35 or */
/*         70 depending on whether N=1 or N .GE. 2.  Under and overflow */
/*         tests are made on the leading term of the asymptotic expansion 
*/
/*         before any extensive computation is done. */

/*     Description of Arguments */

/*         Input */
/*           X      - X .GT. 0.0E0 */
/*           FNU    - order of the initial K function, FNU .GE. 0.0E0 */
/*           KODE   - a parameter to indicate the scaling option */
/*                    KODE=1 returns Y(I)=       K/sub(FNU+I-1)/(X), */
/*                                        I=1,...,N */
/*                    KODE=2 returns Y(I)=EXP(X)*K/sub(FNU+I-1)/(X), */
/*                                        I=1,...,N */
/*           N      - number of members in the sequence, N .GE. 1 */

/*         Output */
/*           y      - a vector whose first n components contain values */
/*                    for the sequence */
/*                    Y(I)=       K/sub(FNU+I-1)/(X), I=1,...,N  or */
/*                    Y(I)=EXP(X)*K/sub(FNU+I-1)/(X), I=1,...,N */
/*                    depending on KODE */
/*           NZ     - number of components of Y set to zero due to */
/*                    underflow with KODE=1, */
/*                    NZ=0   , normal return, computation completed */
/*                    NZ .NE. 0, first NZ components of Y set to zero */
/*                             due to underflow, Y(I)=0.0E0, I=1,...,NZ */

/*     Error Conditions */
/*         Improper input arguments - a fatal error */
/*         Overflow - a fatal error */
/*         Underflow with KODE=1 -  a non-fatal error (NZ .NE. 0) */

/* ***REFERENCES  F. W. J. Olver, Tables of Bessel Functions of Moderate 
*/
/*                 or Large Orders, NPL Mathematical Tables 6, Her */
/*                 Majesty's Stationery Office, London, 1962. */
/*               N. M. Temme, On the numerical evaluation of the modified 
*/
/*                 Bessel function of the third kind, Journal of */
/*                 Computational Physics 19, (1975), pp. 324-337. */
/* ***ROUTINES CALLED  ASYIK, BESK0, BESK0E, BESK1, BESK1E, BESKNU, */
/*                    I1MACH, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790201  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BESK */

    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  BESK */
    /* nn = -i1mach_(&c__12); */
    nn = -i1mach_(&c__15);
    elim = (nn * r1mach_(&c__5) - (float)3.) * (float)2.303;
    xlim = r1mach_(&c__1) * (float)1e3;
    if (*kode < 1 || *kode > 2) {
	goto L280;
    }
    if (*fnu < (float)0.) {
	goto L290;
    }
    if (*x <= (float)0.) {
	goto L300;
    }
    if (*x < xlim) {
	goto L320;
    }
    if (*n < 1) {
	goto L310;
    }
    etx = (real) (*kode - 1);

/*     ND IS A DUMMY VARIABLE FOR N */
/*     GNU IS A DUMMY VARIABLE FOR FNU */
/*     NZ = NUMBER OF UNDERFLOWS ON KODE=1 */

    nd = *n;
    *nz = 0;
    nud = (integer) (*fnu);
    dnu = *fnu - nud;
    gnu = *fnu;
    nn = min(2,nd);
    fn = *fnu + *n - 1;
    fnn = fn;
    if (fn < (float)2.) {
	goto L150;
    }

/*     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION) */
/*     FOR THE LAST ORDER, FNU+N-1.GE.NULIM */

    zn = *x / fn;
    if (zn == (float)0.) {
	goto L320;
    }

    rtz = sqrt(zn * zn + (float)1.);
    gln = log((rtz + (float)1.) / zn);
    t = rtz * ((float)1. - etx) + etx / (zn + rtz);
    cn = -fn * (t - gln);
    if (cn > elim) {
	goto L320;
    }
    if (nud < nulim[nn - 1]) {
	goto L30;
    }
    if (nn == 1) {
	goto L20;
    }

L10:

/*     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION) */
/*     FOR THE FIRST ORDER, FNU.GE.NULIM */

    fn = gnu;
    zn = *x / fn;
    rtz = sqrt(zn * zn + (float)1.);
    gln = log((rtz + (float)1.) / zn);
    t = rtz * ((float)1. - etx) + etx / (zn + rtz);
    cn = -fn * (t - gln);
L20:
    if (cn < -elim) {
	goto L230;
    }

/*     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1.GE.NULIM */

    flgik = (float)-1.;
    asyik_(x, &gnu, kode, &flgik, &rtz, &cn, &nn, &y[1]);
    if (nn == 1) {
	goto L240;
    }
    trx = (float)2. / *x;
    tm = (gnu + gnu + (float)2.) / *x;
    goto L130;

L30:
    if (*kode == 2) {
	goto L40;
    }

/*     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION IN X) 
*/
/*     FOR ORDER DNU */

    if (*x > elim) {
	goto L230;
    }
L40:
    if (dnu != (float)0.) {
	goto L80;
    }
    if (*kode == 2) {
	goto L50;
    }
    s1 = besk0_(x);
    goto L60;
L50:
    s1 = besk0e_(x);
L60:
    if (nud == 0 && nd == 1) {
	goto L120;
    }
    if (*kode == 2) {
	goto L70;
    }
    s2 = besk1_(x);
    goto L90;
L70:
    s2 = besk1e_(x);
    goto L90;
L80:
    nb = 2;
    if (nud == 0 && nd == 1) {
	nb = 1;
    }
    besknu_(x, &dnu, kode, &nb, w, nz);
    s1 = w[0];
    if (nb == 1) {
	goto L120;
    }
    s2 = w[1];
L90:
    trx = (float)2. / *x;
    tm = (dnu + dnu + (float)2.) / *x;
/*     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2) */
    if (nd == 1) {
	--nud;
    }
    if (nud > 0) {
	goto L100;
    }
    if (nd > 1) {
	goto L120;
    }
    s1 = s2;
    goto L120;
L100:
    i__1 = nud;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s = s2;
	s2 = tm * s2 + s1;
	s1 = s;
	tm += trx;
/* L110: */
    }
    if (nd == 1) {
	s1 = s2;
    }
L120:
    y[1] = s1;
    if (nd == 1) {
	goto L240;
    }
    y[2] = s2;
L130:
    if (nd == 2) {
	goto L240;
    }

/*     FORWARD RECUR FROM FNU+2 TO FNU+N-1 */
    i__1 = nd;
    for (i__ = 3; i__ <= i__1; ++i__) {
	y[i__] = tm * y[i__ - 1] + y[i__ - 2];
	tm += trx;
/* L140: */
    }
    goto L240;

L150:
/*     UNDERFLOW TEST FOR KODE=1 */
    if (*kode == 2) {
	goto L160;
    }
    if (*x > elim) {
	goto L230;
    }
L160:
/*     OVERFLOW TEST */
    if (fn <= (float)1.) {
	goto L170;
    }
    if (-fn * (log(*x) - (float).693) > elim) {
	goto L320;
    }
L170:
    if (dnu == (float)0.) {
	goto L180;
    }
    besknu_(x, fnu, kode, &nd, &y[1], &mz);
    goto L240;
L180:
    j = nud;
    if (j == 1) {
	goto L210;
    }
    ++j;
    if (*kode == 2) {
	goto L190;
    }
    y[j] = besk0_(x);
    goto L200;
L190:
    y[j] = besk0e_(x);
L200:
    if (nd == 1) {
	goto L240;
    }
    ++j;
L210:
    if (*kode == 2) {
	goto L220;
    }
    y[j] = besk1_(x);
    goto L240;
L220:
    y[j] = besk1e_(x);
    goto L240;

/*     UPDATE PARAMETERS ON UNDERFLOW */

L230:
    ++nud;
    --nd;
    if (nd == 0) {
	goto L240;
    }
    nn = min(2,nd);
    gnu += (float)1.;
    if (fnn < (float)2.) {
	goto L230;
    }
    if (nud < nulim[nn - 1]) {
	goto L230;
    }
    goto L10;
L240:
    *nz = *n - nd;
    if (*nz == 0) {
	return 0;
    }
    if (nd == 0) {
	goto L260;
    }
    i__1 = nd;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = *n - i__ + 1;
	k = nd - i__ + 1;
	y[j] = y[k];

/* L250: */
    }
L260:
    i__1 = *nz;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = (float)0.;
/* L270: */
    }
    return 0;



L280:
    xermsg_("SLATEC", "BESK", "SCALING OPTION, KODE, NOT 1 OR 2", &c__2, &
	    c__1, 6L, 4L, 32L);
    return 0;
L290:
    xermsg_("SLATEC", "BESK", "ORDER, FNU, LESS THAN ZERO", &c__2, &c__1, 6L, 
	    4L, 26L);
    return 0;
L300:
    xermsg_("SLATEC", "BESK", "X LESS THAN OR EQUAL TO ZERO", &c__2, &c__1, 
	    6L, 4L, 28L);
    return 0;
L310:
    xermsg_("SLATEC", "BESK", "N LESS THAN ONE", &c__2, &c__1, 6L, 4L, 15L);
    return 0;
L320:
    xermsg_("SLATEC", "BESK", "OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL", &
	    c__6, &c__1, 6L, 4L, 43L);
    return -1;
} /* besk_ */

