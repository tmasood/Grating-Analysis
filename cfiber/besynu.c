/* besynu.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__3 = 3;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK BESYNU */
/* Subroutine */ int besynu_(x, fnu, n, y)
real *x, *fnu;
integer *n;
real *y;
{
    /* Initialized data */

    static real x1 = (float)3.;
    static real x2 = (float)20.;
    static real pi = (float)3.14159265358979;
    static real rthpi = (float).797884560802865;
    static real hpi = (float)1.5707963267949;
    static real cc[8] = { (float).577215664901533,(float)-.0420026350340952,(
	    float)-.0421977345555443,(float).007218943246663,(float)
	    -2.152416741149e-4,(float)-2.01348547807e-5,(float)1.133027232e-6,
	    (float)6.116095e-9 };

    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Builtin functions */
    double log(), sin(), sinh(), cosh(), exp(), sqrt(), cos();

    /* Local variables */
    static real coef, relb, flrx, a[120], f, g;
    static integer i__, j, k;
    static real p, q;
    extern doublereal gamma_();
    static real s, a1, a2, etest, g1, g2, s1, s2, t1, t2;
    extern doublereal r1mach_();
    static real cb[120], fc, ak, bk, ck, fk, fn, rb[120];
    static integer kk;
    static real cs, sa, sb, cx;
    static integer nn;
    static real tb, fx, tm, pt, rs, ss, st, rx, cp1, cp2, cs1, cs2;
    extern /* Subroutine */ int xermsg_();
    static real rp1, rp2, rs1, rs2, cbk, cck, arg, rbk, rck, fhs, fks, cpt, 
	    dnu, fmu;
    static integer inu;
    static real tol, etx, smu, rpt, dnu2;

/* ***BEGIN PROLOGUE  BESYNU */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BESY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BESYNU-S, DBSYNU-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         BESYNU computes N member sequences of Y Bessel functions */
/*         Y/SUB(FNU+I-1)/(X), I=1,N for non-negative orders FNU and */
/*         positive X. Equations of the references are implemented on */
/*         small orders DNU for Y/SUB(DNU)/(X) and Y/SUB(DNU+1)/(X). */
/*         Forward recursion with the three term recursion relation */
/*         generates higher orders FNU+I-1, I=1,...,N. */

/*         To start the recursion FNU is normalized to the interval */
/*         -0.5.LE.DNU.LT.0.5. A special form of the power series is */
/*         implemented on 0.LT.X.LE.X1 while the Miller algorithm for the 
*/
/*         K Bessel function in terms of the confluent hypergeometric */
/*         function U(FNU+0.5,2*FNU+1,I*X) is implemented on X1.LT.X.LE.X 
*/
/*         Here I is the complex number SQRT(-1.). */
/*         For X.GT.X2, the asymptotic expansion for large X is used. */
/*         When FNU is a half odd integer, a special formula for */
/*         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion. */

/*         BESYNU assumes that a significant digit SINH(X) function is */
/*         available. */

/*     Description of Arguments */

/*         Input */
/*           X      - X.GT.0.0E0 */
/*           FNU    - Order of initial Y function, FNU.GE.0.0E0 */
/*           N      - Number of members of the sequence, N.GE.1 */

/*         Output */
/*           Y      - A vector whose first N components contain values */
/*                    for the sequence Y(I)=Y/SUB(FNU+I-1), I=1,N. */

/*     Error Conditions */
/*         Improper input arguments - a fatal error */
/*         Overflow - a fatal error */

/* ***SEE ALSO  BESY */
/* ***REFERENCES  N. M. Temme, On the numerical evaluation of the ordinary
 */
/*                 Bessel function of the second kind, Journal of */
/*                 Computational Physics 21, (1976), pp. 343-350. */
/*               N. M. Temme, On the numerical evaluation of the modified 
*/
/*                 Bessel function of the third kind, Journal of */
/*                 Computational Physics 19, (1975), pp. 324-337. */
/* ***ROUTINES CALLED  GAMMA, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/*   910408  Updated the AUTHOR and REFERENCES sections.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BESYNU */

    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  BESYNU */
    ak = r1mach_(&c__3);
    tol = dmax(ak,(float)1e-15);
    if (*x <= (float)0.) {
	goto L270;
    }
    if (*fnu < (float)0.) {
	goto L280;
    }
    if (*n < 1) {
	goto L290;
    }
    rx = (float)2. / *x;
    inu = (integer) (*fnu + (float).5);
    dnu = *fnu - inu;
    if (dabs(dnu) == (float).5) {
	goto L260;
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
    g1 = -(s + s);
    goto L50;
L40:
    g1 = (t1 - t2) / dnu;
L50:
    g2 = t1 + t2;
    smu = (float)1.;
    fc = (float)1. / pi;
    flrx = log(rx);
    fmu = dnu * flrx;
    tm = (float)0.;
    if (dnu == (float)0.) {
	goto L60;
    }
    tm = sin(dnu * hpi) / dnu;
    tm = (dnu + dnu) * tm * tm;
    fc = dnu / sin(dnu * pi);
    if (fmu != (float)0.) {
	smu = sinh(fmu) / fmu;
    }
L60:
    f = fc * (g1 * cosh(fmu) + g2 * flrx * smu);
    fx = exp(fmu);
    p = fc * t1 * fx;
    q = fc * t2 / fx;
    g = f + tm * q;
    ak = (float)1.;
    ck = (float)1.;
    bk = (float)1.;
    s1 = g;
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
    g = f + tm * q;
    ck = -ck * cx / ak;
    t1 = ck * g;
    s1 += t1;
    bk = bk + ak + ak + (float)1.;
    ak += (float)1.;
    s = dabs(t1) / (dabs(s1) + (float)1.);
    if (s > tol) {
	goto L70;
    }
L80:
    y[1] = -s1;
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
    g = f + tm * q;
    ck = -ck * cx / ak;
    t1 = ck * g;
    s1 += t1;
    t2 = ck * (p - ak * g);
    s2 += t2;
    bk = bk + ak + ak + (float)1.;
    ak += (float)1.;
    s = dabs(t1) / (dabs(s1) + (float)1.) + dabs(t2) / (dabs(s2) + (float)1.);
    if (s > tol) {
	goto L100;
    }
L110:
    s2 = -s2 * rx;
    s1 = -s1;
    goto L160;
L120:
    coef = rthpi / sqrt(*x);
    if (*x > x2) {
	goto L210;
    }

/*     MILLER ALGORITHM FOR X1.LT.X.LE.X2 */

    etest = cos(pi * dnu) / (pi * *x * tol);
    fks = (float)1.;
    fhs = (float).25;
    fk = (float)0.;
    rck = (float)2.;
    cck = *x + *x;
    rp1 = (float)0.;
    cp1 = (float)0.;
    rp2 = (float)1.;
    cp2 = (float)0.;
    k = 0;
L130:
    ++k;
    fk += (float)1.;
    ak = (fhs - dnu2) / (fks + fk);
    pt = fk + (float)1.;
    rbk = rck / pt;
    cbk = cck / pt;
    rpt = rp2;
    cpt = cp2;
    rp2 = rbk * rpt - cbk * cpt - ak * rp1;
    cp2 = cbk * rpt + rbk * cpt - ak * cp1;
    rp1 = rpt;
    cp1 = cpt;
    rb[k - 1] = rbk;
    cb[k - 1] = cbk;
    a[k - 1] = ak;
    rck += (float)2.;
    fks = fks + fk + fk + (float)1.;
    fhs = fhs + fk + fk;
/* Computing MAX */
    r__1 = dabs(rp1), r__2 = dabs(cp1);
    pt = dmax(r__1,r__2);
/* Computing 2nd power */
    r__1 = rp1 / pt;
/* Computing 2nd power */
    r__2 = cp1 / pt;
    fc = r__1 * r__1 + r__2 * r__2;
    pt = pt * sqrt(fc) * fk;
    if (etest > pt) {
	goto L130;
    }
    kk = k;
    rs = (float)1.;
    cs = (float)0.;
    rp1 = (float)0.;
    cp1 = (float)0.;
    rp2 = (float)1.;
    cp2 = (float)0.;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rpt = rp2;
	cpt = cp2;
	rp2 = (rb[kk - 1] * rpt - cb[kk - 1] * cpt - rp1) / a[kk - 1];
	cp2 = (cb[kk - 1] * rpt + rb[kk - 1] * cpt - cp1) / a[kk - 1];
	rp1 = rpt;
	cp1 = cpt;
	rs += rp2;
	cs += cp2;
	--kk;
/* L140: */
    }
/* Computing MAX */
    r__1 = dabs(rs), r__2 = dabs(cs);
    pt = dmax(r__1,r__2);
/* Computing 2nd power */
    r__1 = rs / pt;
/* Computing 2nd power */
    r__2 = cs / pt;
    fc = r__1 * r__1 + r__2 * r__2;
    pt *= sqrt(fc);
    rs1 = (rp2 * (rs / pt) + cp2 * (cs / pt)) / pt;
    cs1 = (cp2 * (rs / pt) - rp2 * (cs / pt)) / pt;
    fc = hpi * (dnu - (float).5) - *x;
    p = cos(fc);
    q = sin(fc);
    s1 = (cs1 * q - rs1 * p) * coef;
    if (inu > 0 || *n > 1) {
	goto L150;
    }
    y[1] = s1;
    return 0;
L150:
/* Computing MAX */
    r__1 = dabs(rp2), r__2 = dabs(cp2);
    pt = dmax(r__1,r__2);
/* Computing 2nd power */
    r__1 = rp2 / pt;
/* Computing 2nd power */
    r__2 = cp2 / pt;
    fc = r__1 * r__1 + r__2 * r__2;
    pt *= sqrt(fc);
    rpt = dnu + (float).5 - (rp1 * (rp2 / pt) + cp1 * (cp2 / pt)) / pt;
    cpt = *x - (cp1 * (rp2 / pt) - rp1 * (cp2 / pt)) / pt;
    cs2 = cs1 * cpt - rs1 * rpt;
    rs2 = rpt * cs1 + rs1 * cpt;
    s2 = (rs2 * q + cs2 * p) * coef / *x;

/*     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION */

L160:
    ck = (dnu + dnu + (float)2.) / *x;
    if (*n == 1) {
	--inu;
    }
    if (inu > 0) {
	goto L170;
    }
    if (*n > 1) {
	goto L190;
    }
    s1 = s2;
    goto L190;
L170:
    i__1 = inu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	st = s2;
	s2 = ck * s2 - s1;
	s1 = st;
	ck += rx;
/* L180: */
    }
    if (*n == 1) {
	s1 = s2;
    }
L190:
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
	y[i__] = ck * y[i__ - 1] - y[i__ - 2];
	ck += rx;
/* L200: */
    }
    return 0;

/*     ASYMPTOTIC EXPANSION FOR LARGE X, X.GT.X2 */

L210:
    nn = 2;
    if (inu == 0 && *n == 1) {
	nn = 1;
    }
    dnu2 = dnu + dnu;
    fmu = (float)0.;
    if (dabs(dnu2) < tol) {
	goto L220;
    }
    fmu = dnu2 * dnu2;
L220:
    arg = *x - hpi * (dnu + (float).5);
    sa = sin(arg);
    sb = cos(arg);
    etx = *x * (float)8.;
    i__1 = nn;
    for (k = 1; k <= i__1; ++k) {
	s1 = s2;
	t2 = (fmu - (float)1.) / etx;
	ss = t2;
	relb = tol * dabs(t2);
	t1 = etx;
	s = (float)1.;
	fn = (float)1.;
	ak = (float)0.;
	for (j = 1; j <= 13; ++j) {
	    t1 += etx;
	    ak += (float)8.;
	    fn += ak;
	    t2 = -t2 * (fmu - fn) / t1;
	    s += t2;
	    t1 += etx;
	    ak += (float)8.;
	    fn += ak;
	    t2 = t2 * (fmu - fn) / t1;
	    ss += t2;
	    if (dabs(t2) <= relb) {
		goto L240;
	    }
/* L230: */
	}
L240:
	s2 = coef * (s * sa + ss * sb);
	fmu = fmu + dnu * (float)8. + (float)4.;
	tb = sa;
	sa = -sb;
	sb = tb;
/* L250: */
    }
    if (nn > 1) {
	goto L160;
    }
    s1 = s2;
    goto L190;

/*     FNU=HALF ODD INTEGER CASE */

L260:
    coef = rthpi / sqrt(*x);
    s1 = coef * sin(*x);
    s2 = -coef * cos(*x);
    goto L160;


L270:
    xermsg_("SLATEC", "BESYNU", "X NOT GREATER THAN ZERO", &c__2, &c__1, 6L, 
	    6L, 23L);
    return 0;
L280:
    xermsg_("SLATEC", "BESYNU", "FNU NOT ZERO OR POSITIVE", &c__2, &c__1, 6L, 
	    6L, 24L);
    return 0;
L290:
    xermsg_("SLATEC", "BESYNU", "N NOT GREATER THAN 0", &c__2, &c__1, 6L, 6L, 
	    20L);
    return 0;
} /* besynu_ */

