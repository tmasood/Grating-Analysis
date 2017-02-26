/* asyik.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__3 = 3;

/* DECK ASYIK */
/* Subroutine */ int asyik_(x, fnu, kode, flgik, ra, arg, in, y)
real *x, *fnu;
integer *kode;
real *flgik, *ra, *arg;
integer *in;
real *y;
{
    /* Initialized data */

    static real con[2] = { (float).398942280401432678,(float)
	    1.25331413731550025 };
    static real c__[65] = { (float)-.208333333333333,(float).125,(float)
	    .334201388888889,(float)-.401041666666667,(float).0703125,(float)
	    -1.02581259645062,(float)1.84646267361111,(float)-.8912109375,(
	    float).0732421875,(float)4.66958442342625,(float)-11.207002616223,
	    (float)8.78912353515625,(float)-2.3640869140625,(float)
	    .112152099609375,(float)-28.2120725582002,(float)84.6362176746007,
	    (float)-91.81824154324,(float)42.5349987453885,(float)
	    -7.36879435947963,(float).227108001708984,(float)212.570130039217,
	    (float)-765.252468141182,(float)1059.990452528,(float)
	    -699.579627376133,(float)218.190511744212,(float)
	    -26.4914304869516,(float).572501420974731,(float)
	    -1919.45766231841,(float)8061.72218173731,(float)
	    -13586.5500064341,(float)11655.3933368645,(float)-5305.6469786134,
	    (float)1200.90291321635,(float)-108.090919788395,(float)
	    1.72772750258446,(float)20204.2913309661,(float)-96980.5983886375,
	    (float)192547.001232532,(float)-203400.177280416,(float)
	    122200.464983017,(float)-41192.6549688976,(float)7109.51430248936,
	    (float)-493.915304773088,(float)6.07404200127348,(float)
	    -242919.187900551,(float)1311763.61466298,(float)
	    -2998015.91853811,(float)3763271.2976564,(float)-2813563.22658653,
	    (float)1268365.27332162,(float)-331645.172484564,(float)
	    45218.7689813627,(float)-2499.83048181121,(float)24.3805296995561,
	    (float)3284469.85307204,(float)-19706819.1184322,(float)
	    50952602.4926646,(float)-74105148.2115327,(float)66344512.274729,(
	    float)-37567176.6607634,(float)13288767.1664218,(float)
	    -2785618.12808645,(float)308186.404612662,(float)-13886.089753717,
	    (float)110.017140269247 };

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Builtin functions */
    double sqrt(), log(), exp(), r_sign();

    /* Local variables */
    static real coef;
    static integer j, k, l;
    static real t, z__, s1, s2, t2;
    extern doublereal r1mach_();
    static real ak, ap, fn;
    static integer kk, jn;
    static real gln, tol, etx;

/* ***BEGIN PROLOGUE  ASYIK */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BESI and BESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (ASYIK-S, DASYIK-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*                    ASYIK computes Bessel functions I and K */
/*                  for arguments X.GT.0.0 and orders FNU.GE.35 */
/*                  on FLGIK = 1 and FLGIK = -1 respectively. */

/*                                    INPUT */

/*      X    - argument, X.GT.0.0E0 */
/*      FNU  - order of first Bessel function */
/*      KODE - a parameter to indicate the scaling option */
/*             KODE=1 returns Y(I)=        I/SUB(FNU+I-1)/(X), I=1,IN */
/*                    or      Y(I)=        K/SUB(FNU+I-1)/(X), I=1,IN */
/*                    on FLGIK = 1.0E0 or FLGIK = -1.0E0 */
/*             KODE=2 returns Y(I)=EXP(-X)*I/SUB(FNU+I-1)/(X), I=1,IN */
/*                    or      Y(I)=EXP( X)*K/SUB(FNU+I-1)/(X), I=1,IN */
/*                    on FLGIK = 1.0E0 or FLGIK = -1.0E0 */
/*     FLGIK - selection parameter for I or K function */
/*             FLGIK =  1.0E0 gives the I function */
/*             FLGIK = -1.0E0 gives the K function */
/*        RA - SQRT(1.+Z*Z), Z=X/FNU */
/*       ARG - argument of the leading exponential */
/*        IN - number of functions desired, IN=1 or 2 */

/*                                    OUTPUT */

/*         Y - a vector whose first in components contain the sequence */

/*     Abstract */
/*         ASYIK implements the uniform asymptotic expansion of */
/*         the I and K Bessel functions for FNU.GE.35 and real */
/*         X.GT.0.0E0. The forms are identical except for a change */
/*         in sign of some of the terms. This change in sign is */
/*         accomplished by means of the flag FLGIK = 1 or -1. */

/* ***SEE ALSO  BESI, BESK */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR section.  (WRB) */
/* ***END PROLOGUE  ASYIK */

    /* Parameter adjustments */
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ASYIK */
    tol = r1mach_(&c__3);
    tol = dmax(tol,(float)1e-15);
    fn = *fnu;
    z__ = ((float)3. - *flgik) / (float)2.;
    kk = (integer) z__;
    i__1 = *in;
    for (jn = 1; jn <= i__1; ++jn) {
	if (jn == 1) {
	    goto L10;
	}
	fn -= *flgik;
	z__ = *x / fn;
	*ra = sqrt(z__ * z__ + (float)1.);
	gln = log((*ra + (float)1.) / z__);
	etx = (real) (*kode - 1);
	t = *ra * ((float)1. - etx) + etx / (z__ + *ra);
	*arg = fn * (t - gln) * *flgik;
L10:
	coef = exp(*arg);
	t = (float)1. / *ra;
	t2 = t * t;
	t /= fn;
	t = r_sign(&t, flgik);
	s2 = (float)1.;
	ap = (float)1.;
	l = 0;
	for (k = 2; k <= 11; ++k) {
	    ++l;
	    s1 = c__[l - 1];
	    i__2 = k;
	    for (j = 2; j <= i__2; ++j) {
		++l;
		s1 = s1 * t2 + c__[l - 1];
/* L20: */
	    }
	    ap *= t;
	    ak = ap * s1;
	    s2 += ak;
/* Computing MAX */
	    r__1 = dabs(ak), r__2 = dabs(ap);
	    if (dmax(r__1,r__2) < tol) {
		goto L40;
	    }
/* L30: */
	}
L40:
	t = dabs(t);
	y[jn] = s2 * coef * sqrt(t) * con[kk - 1];
/* L50: */
    }
    return 0;
} /* asyik_ */

