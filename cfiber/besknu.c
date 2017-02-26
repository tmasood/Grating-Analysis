/* besknu.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__12 = 12;
static integer c__5 = 5;
static integer c__3 = 3;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK BESKNU */
/* Subroutine */ int besknu_(x, fnu, kode, n, y, nz)
real *x, *fnu;
integer *kode, *n;
real *y;
integer *nz;
{
    /* Initialized data */

    static real x1 = (float)2.;
    static real x2 = (float)17.;
    static real pi = (float)3.14159265358979;
    static real rthpi = (float)1.2533141373155;
    static real cc[8] = { (float).577215664901533,(float)-.0420026350340952,(
	    float)-.0421977345555443,(float).007218943246663,(float)
	    -2.152416741149e-4,(float)-2.01348547807e-5,(float)1.133027232e-6,
	    (float)6.116095e-9 };

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(), sin(), sinh(), cosh(), exp(), sqrt(), cos();

    /* Local variables */
    static real coef, elim, flrx, a[160], b[160], f;
    static integer i__, j, k;
    static real p, q;
    static integer iflag;
    static real s;
    extern doublereal gamma_();
    static integer koded;
    static real a1, a2, etest, g1, g2, p1;
    extern integer i1mach_();
    static real p2, s1, s2, t1, t2;
    extern doublereal r1mach_();
    static real fc, ak, bk, ck, dk, fk;
    static integer kk;
    static real cx;
    static integer nn;
    static real ex, tm, pt, st, rx;
    extern /* Subroutine */ int xermsg_();
    static real fhs, fks, dnu, fmu;
    static integer inu;
    static real sqk, tol, smu, dnu2;

/* ***BEGIN PROLOGUE  BESKNU */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BESKNU-S, DBSKNU-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         BESKNU computes N member sequences of K Bessel functions */
/*         K/SUB(FNU+I-1)/(X), I=1,N for non-negative orders FNU and */
/*         positive X. Equations of the references are implemented on */
/*         small orders DNU for K/SUB(DNU)/(X) and K/SUB(DNU+1)/(X). */
/*         Forward recursion with the three term recursion relation */
/*         generates higher orders FNU+I-1, I=1,...,N. The parameter */
/*         KODE permits K/SUB(FNU+I-1)/(X) values or scaled values */
/*         EXP(X)*K/SUB(FNU+I-1)/(X), I=1,N to be returned. */

/*         To start the recursion FNU is normalized to the interval */
/*         -0.5.LE.DNU.LT.0.5. A special form of the power series is */
/*         implemented on 0.LT.X.LE.X1 while the Miller algorithm for the 
*/
/*         K Bessel function in terms of the confluent hypergeometric */
/*         function U(FNU+0.5,2*FNU+1,X) is implemented on X1.LT.X.LE.X2. 
*/
/*         For X.GT.X2, the asymptotic expansion for large X is used. */
/*         When FNU is a half odd integer, a special formula for */
/*         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion. */

/*         BESKNU assumes that a significant digit SINH(X) function is */
/*         available. */

/*     Description of Arguments */

/*         Input */
/*           X      - X.GT.0.0E0 */
/*           FNU    - Order of initial K function, FNU.GE.0.0E0 */
/*           N      - Number of members of the sequence, N.GE.1 */
/*           KODE   - A parameter to indicate the scaling option */
/*                    KODE= 1  returns */
/*                             Y(I)=       K/SUB(FNU+I-1)/(X) */
/*                                  I=1,...,N */
/*                        = 2  returns */
/*                             Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X) */
/*                                  I=1,...,N */

/*         Output */
/*           Y      - A vector whose first N components contain values */
/*                    for the sequence */
/*                    Y(I)=       K/SUB(FNU+I-1)/(X), I=1,...,N or */
/*                    Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N */
/*                    depending on KODE */
/*           NZ     - Number of components set to zero due to */
/*                    underflow, */
/*                    NZ= 0   , Normal return */
/*                    NZ.NE.0 , First NZ components of Y set to zero */
/*                              due to underflow, Y(I)=0.0E0,I=1,...,NZ */

/*     Error Conditions */
/*         Improper input arguments - a fatal error */
/*         Overflow - a fatal error */
/*         Underflow with KODE=1 - a non-fatal error (NZ.NE.0) */

/* ***SEE ALSO  BESK */
/* ***REFERENCES  N. M. Temme, On the numerical evaluation of the modified
 */
/*                 Bessel function of the third kind, Journal of */
/*                 Computational Physics 19, (1975), pp. 324-337. */
/* ***ROUTINES CALLED  GAMMA, I1MACH, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790201  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/*   910408  Updated the AUTHOR and REFERENCES sections.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BESKNU */

    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  BESKNU */
    kk = -i1mach_(&c__12);
    elim = (kk * r1mach_(&c__5) - (float)3.) * (float)2.303;
    ak = r1mach_(&c__3);
    tol = dmax(ak,(float)1e-15);
    if (*x <= (float)0.) {
	goto L350;
    }
    if (*fnu < (float)0.) {
	goto L360;
    }
    if (*kode < 1 || *kode > 2) {
	goto L370;
    }
    if (*n < 1) {
	goto L380;
    }
    *nz = 0;
    iflag = 0;
    koded = *kode;
    rx = (float)2. / *x;
    inu = (integer) (*fnu + (float).5);
    dnu = *fnu - inu;
    if (dabs(dnu) == (float).5) {
	goto L120;
    }
    dnu2 = (float)0.;
    if (dabs(dnu) < tol) {
	goto L10;
    }
    dnu2 = dnu * dnu;
L10:
    if (*x > x1) {
	goto L120;
    }

/*     SERIES FOR X.LE.X1 */

    a1 = (float)1. - dnu;
    a2 = dnu + (float)1.;
    t1 = (float)1. / gamma_(&a1);
    t2 = (float)1. / gamma_(&a2);
    if (dabs(dnu) > (float).1) {
	goto L40;
    }
/*     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU) */
    s = cc[0];
    ak = (float)1.;
    for (k = 2; k <= 8; ++k) {
	ak *= dnu2;
	tm = cc[k - 1] * ak;
	s += tm;
	if (dabs(tm) < tol) {
	    goto L30;
	}
/* L20: */
    }
L30:
    g1 = -s;
    goto L50;
L40:
    g1 = (t1 - t2) / (dnu + dnu);
L50:
    g2 = (t1 + t2) * (float).5;
    smu = (float)1.;
    fc = (float)1.;
    flrx = log(rx);
    fmu = dnu * flrx;
    if (dnu == (float)0.) {
	goto L60;
    }
    fc = dnu * pi;
    fc /= sin(fc);
    if (fmu != (float)0.) {
	smu = sinh(fmu) / fmu;
    }
L60:
    f = fc * (g1 * cosh(fmu) + g2 * flrx * smu);
    fc = exp(fmu);
    p = fc * (float).5 / t2;
    q = (float).5 / (fc * t1);
    ak = (float)1.;
    ck = (float)1.;
    bk = (float)1.;
    s1 = f;
    s2 = p;
    if (inu > 0 || *n > 1) {
	goto L90;
    }
    if (*x < tol) {
	goto L80;
    }
    cx = *x * *x * (float).25;
L70:
    f = (ak * f + p + q) / (bk - dnu2);
    p /= ak - dnu;
    q /= ak + dnu;
    ck = ck * cx / ak;
    t1 = ck * f;
    s1 += t1;
    bk = bk + ak + ak + (float)1.;
    ak += (float)1.;
    s = dabs(t1) / (dabs(s1) + (float)1.);
    if (s > tol) {
	goto L70;
    }
L80:
    y[1] = s1;
    if (koded == 1) {
	return 0;
    }
    y[1] = s1 * exp(*x);
    return 0;
L90:
    if (*x < tol) {
	goto L110;
    }
    cx = *x * *x * (float).25;
L100:
    f = (ak * f + p + q) / (bk - dnu2);
    p /= ak - dnu;
    q /= ak + dnu;
    ck = ck * cx / ak;
    t1 = ck * f;
    s1 += t1;
    t2 = ck * (p - ak * f);
    s2 += t2;
    bk = bk + ak + ak + (float)1.;
    ak += (float)1.;
    s = dabs(t1) / (dabs(s1) + (float)1.) + dabs(t2) / (dabs(s2) + (float)1.);
    if (s > tol) {
	goto L100;
    }
L110:
    s2 *= rx;
    if (koded == 1) {
	goto L170;
    }
    f = exp(*x);
    s1 *= f;
    s2 *= f;
    goto L170;
L120:
    coef = rthpi / sqrt(*x);
    if (koded == 2) {
	goto L130;
    }
    if (*x > elim) {
	goto L330;
    }
    coef *= exp(-(*x));
L130:
    if (dabs(dnu) == (float).5) {
	goto L340;
    }
    if (*x > x2) {
	goto L280;
    }

/*     MILLER ALGORITHM FOR X1.LT.X.LE.X2 */

    etest = cos(pi * dnu) / (pi * *x * tol);
    fks = (float)1.;
    fhs = (float).25;
    fk = (float)0.;
    ck = *x + *x + (float)2.;
    p1 = (float)0.;
    p2 = (float)1.;
    k = 0;
L140:
    ++k;
    fk += (float)1.;
    ak = (fhs - dnu2) / (fks + fk);
    bk = ck / (fk + (float)1.);
    pt = p2;
    p2 = bk * p2 - ak * p1;
    p1 = pt;
    a[k - 1] = ak;
    b[k - 1] = bk;
    ck += (float)2.;
    fks = fks + fk + fk + (float)1.;
    fhs = fhs + fk + fk;
    if (etest > fk * p1) {
	goto L140;
    }
    kk = k;
    s = (float)1.;
    p1 = (float)0.;
    p2 = (float)1.;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pt = p2;
	p2 = (b[kk - 1] * p2 - p1) / a[kk - 1];
	p1 = pt;
	s += p2;
	--kk;
/* L150: */
    }
    s1 = coef * (p2 / s);
    if (inu > 0 || *n > 1) {
	goto L160;
    }
    goto L200;
L160:
    s2 = s1 * (*x + dnu + (float).5 - p1 / p2) / *x;

/*     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION */

L170:
    ck = (dnu + dnu + (float)2.) / *x;
    if (*n == 1) {
	--inu;
    }
    if (inu > 0) {
	goto L180;
    }
    if (*n > 1) {
	goto L200;
    }
    s1 = s2;
    goto L200;
L180:
    i__1 = inu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	st = s2;
	s2 = ck * s2 + s1;
	s1 = st;
	ck += rx;
/* L190: */
    }
    if (*n == 1) {
	s1 = s2;
    }
L200:
    if (iflag == 1) {
	goto L220;
    }
    y[1] = s1;
    if (*n == 1) {
	return 0;
    }
    y[2] = s2;
    if (*n == 2) {
	return 0;
    }
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; ++i__) {
	y[i__] = ck * y[i__ - 1] + y[i__ - 2];
	ck += rx;
/* L210: */
    }
    return 0;
/*     IFLAG=1 CASES */
L220:
    s = -(*x) + log(s1);
    y[1] = (float)0.;
    *nz = 1;
    if (s < -elim) {
	goto L230;
    }
    y[1] = exp(s);
    *nz = 0;
L230:
    if (*n == 1) {
	return 0;
    }
    s = -(*x) + log(s2);
    y[2] = (float)0.;
    ++(*nz);
    if (s < -elim) {
	goto L240;
    }
    --(*nz);
    y[2] = exp(s);
L240:
    if (*n == 2) {
	return 0;
    }
    kk = 2;
    if (*nz < 2) {
	goto L260;
    }
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; ++i__) {
	kk = i__;
	st = s2;
	s2 = ck * s2 + s1;
	s1 = st;
	ck += rx;
	s = -(*x) + log(s2);
	++(*nz);
	y[i__] = (float)0.;
	if (s < -elim) {
	    goto L250;
	}
	y[i__] = exp(s);
	--(*nz);
	goto L260;
L250:
	;
    }
    return 0;
L260:
    if (kk == *n) {
	return 0;
    }
    s2 = s2 * ck + s1;
    ck += rx;
    ++kk;
    y[kk] = exp(-(*x) + log(s2));
    if (kk == *n) {
	return 0;
    }
    ++kk;
    i__1 = *n;
    for (i__ = kk; i__ <= i__1; ++i__) {
	y[i__] = ck * y[i__ - 1] + y[i__ - 2];
	ck += rx;
/* L270: */
    }
    return 0;

/*     ASYMPTOTIC EXPANSION FOR LARGE X, X.GT.X2 */

/*     IFLAG=0 MEANS NO UNDERFLOW OCCURRED */
/*     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH */
/*     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD */
/*     RECURSION */
L280:
    nn = 2;
    if (inu == 0 && *n == 1) {
	nn = 1;
    }
    dnu2 = dnu + dnu;
    fmu = (float)0.;
    if (dabs(dnu2) < tol) {
	goto L290;
    }
    fmu = dnu2 * dnu2;
L290:
    ex = *x * (float)8.;
    s2 = (float)0.;
    i__1 = nn;
    for (k = 1; k <= i__1; ++k) {
	s1 = s2;
	s = (float)1.;
	ak = (float)0.;
	ck = (float)1.;
	sqk = (float)1.;
	dk = ex;
	for (j = 1; j <= 30; ++j) {
	    ck = ck * (fmu - sqk) / dk;
	    s += ck;
	    dk += ex;
	    ak += (float)8.;
	    sqk += ak;
	    if (dabs(ck) < tol) {
		goto L310;
	    }
/* L300: */
	}
L310:
	s2 = s * coef;
	fmu = fmu + dnu * (float)8. + (float)4.;
/* L320: */
    }
    if (nn > 1) {
	goto L170;
    }
    s1 = s2;
    goto L200;
L330:
    koded = 2;
    iflag = 1;
    goto L120;

/*     FNU=HALF ODD INTEGER CASE */

L340:
    s1 = coef;
    s2 = coef;
    goto L170;


L350:
    xermsg_("SLATEC", "BESKNU", "X NOT GREATER THAN ZERO", &c__2, &c__1, 6L, 
	    6L, 23L);
    return 0;
L360:
    xermsg_("SLATEC", "BESKNU", "FNU NOT ZERO OR POSITIVE", &c__2, &c__1, 6L, 
	    6L, 24L);
    return 0;
L370:
    xermsg_("SLATEC", "BESKNU", "KODE NOT 1 OR 2", &c__2, &c__1, 6L, 6L, 15L);
    return 0;
L380:
    xermsg_("SLATEC", "BESKNU", "N NOT GREATER THAN 0", &c__2, &c__1, 6L, 6L, 
	    20L);
    return 0;
} /* besknu_ */

