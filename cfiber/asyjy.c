/* asyjy.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__3 = 3;
static integer c__5 = 5;
static integer c__12 = 12;
static integer c__11 = 11;
static integer c__1 = 1;

/* DECK ASYJY */
/* Subroutine */ int asyjy_(funjy, x, fnu, flgjy, in, y, wk, iflw)
/* Subroutine */ int (*funjy) ();
real *x, *fnu, *flgjy;
integer *in;
real *y, *wk;
integer *iflw;
{
    /* Initialized data */

    static real con2 = (float).333333333333333;
    static real con548 = (float).104166666666667;
    static real ar[8] = { (float).0835503472222222,(float).128226574556327,(
	    float).29184902646414,(float).881627267443758,(float)
	    3.32140828186277,(float)14.9957629868626,(float)78.9230130115865,(
	    float)474.451538868264 };
    static real br[10] = { (float)-.145833333333333,(float)-.0987413194444444,
	    (float)-.143312053915895,(float)-.317227202678414,(float)
	    -.94242914795712,(float)-3.51120304082635,(float)-15.727263620368,
	    (float)-82.2814390971859,(float)-492.355370523671,(float)
	    -3316.21856854797 };
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
    static real gama[26] = { (float).629960524947437,(float).251984209978975,(
	    float).154790300415656,(float).110713062416159,(float)
	    .0857309395527395,(float).0697161316958684,(float)
	    .0586085671893714,(float).0504698873536311,(float)
	    .0442600580689155,(float).039372066154351,(float)
	    .0354283195924455,(float).0321818857502098,(float)
	    .0294646240791158,(float).0271581677112934,(float)
	    .0251768272973862,(float).0234570755306079,(float)
	    .0219508390134907,(float).0206210828235646,(float)
	    .0194388240897881,(float).0183810633800683,(float)
	    .0174293213231963,(float).0165685837786612,(float)
	    .0157865285987918,(float).0150729501494096,(float)
	    .0144193250839955,(float).0138184805735342 };
    static real tols = (float)-6.90775527898214;
    static real con1 = (float).666666666666667;
    static struct {
	real e_1[104];
	} equiv_1 = { (float)-.00444444444444444, (float)-9.22077922077922e-4,
		 (float)-8.84892884892885e-5, (float)1.6592768783245e-4, (
		float)2.46691372741793e-4, (float)2.65995589346255e-4, (float)
		2.61824297061501e-4, (float)2.48730437344656e-4, (float)
		2.32721040083232e-4, (float)2.16362485712365e-4, (float)
		2.00738858762752e-4, (float)1.86267636637545e-4, (float)
		1.73060775917876e-4, (float)1.61091705929016e-4, (float)
		1.50274774160908e-4, (float)1.4050349739127e-4, (float)
		1.31668816545923e-4, (float)1.23667445598253e-4, (float)
		1.16405271474738e-4, (float)1.09798298372713e-4, (float)
		1.03772410422993e-4, (float)9.82626078369363e-5, (float)
		9.32120517249503e-5, (float)8.85710852478712e-5, (float)
		8.429631057157e-5, (float)8.03497548407791e-5, (float)
		6.93735541354589e-4, (float)2.32241745182922e-4, (float)
		-1.41986273556691e-5, (float)-1.16444931672049e-4, (float)
		-1.50803558053049e-4, (float)-1.55121924918096e-4, (float)
		-1.46809756646466e-4, (float)-1.33815503867491e-4, (float)
		-1.19744975684254e-4, (float)-1.06184319207974e-4, (float)
		-9.37699549891194e-5, (float)-8.26923045588193e-5, (float)
		-7.29374348155221e-5, (float)-6.44042357721016e-5, (float)
		-5.69611566009369e-5, (float)-5.04731044303562e-5, (float)
		-4.48134868008883e-5, (float)-3.98688727717599e-5, (float)
		-3.55400532972042e-5, (float)-3.17414256609022e-5, (float)
		-2.83996793904175e-5, (float)-2.54522720634871e-5, (float)
		-2.28459297164725e-5, (float)-2.05352753106481e-5, (float)
		-1.84816217627666e-5, (float)-1.66519330021394e-5, (float)
		-3.54211971457744e-4, (float)-1.56161263945159e-4, (float)
		3.04465503594936e-5, (float)1.30198655773243e-4, (float)
		1.67471106699712e-4, (float)1.70222587683593e-4, (float)
		1.56501427608595e-4, (float)1.36339170977445e-4, (float)
		1.14886692029825e-4, (float)9.45869093034688e-5, (float)
		7.64498419250898e-5, (float)6.07570334965197e-5, (float)
		4.74394299290509e-5, (float)3.62757512005344e-5, (float)
		2.69939714979225e-5, (float)1.93210938247939e-5, (float)
		1.30056674793963e-5, (float)7.82620866744497e-6, (float)
		3.59257485819352e-6, (float)1.44040049814252e-7, (float)
		-2.65396769697939e-6, (float)-4.91346867098486e-6, (float)
		-6.72739296091248e-6, (float)-8.17269379678658e-6, (float)
		-9.31304715093561e-6, (float)-1.02011418798016e-5, (float)
		3.78194199201773e-4, (float)2.02471952761816e-4, (float)
		-6.37938506318862e-5, (float)-2.38598230603006e-4, (float)
		-3.10916256027362e-4, (float)-3.13680115247576e-4, (float)
		-2.78950273791323e-4, (float)-2.28564082619141e-4, (float)
		-1.75245280340847e-4, (float)-1.2554406306069e-4, (float)
		-8.22982872820208e-5, (float)-4.62860730588116e-5, (float)
		-1.72334302366962e-5, (float)5.60690482304602e-6, (float)
		2.31395443148287e-5, (float)3.62642745856794e-5, (float)
		4.58006124490189e-5, (float)5.24595294959114e-5, (float)
		5.68396208545815e-5, (float)5.94349820393104e-5, (float)
		6.06478527578422e-5, (float)6.08023907788436e-5, (float)
		6.0157789453946e-5, (float)5.89199657344698e-5, (float)
		5.72515823777593e-5, (float)5.52804375585853e-5 };

    static struct {
	real e_1[130];
	} equiv_4 = { (float).0179988721413553, (float).00559964911064388, (
		float).00288501402231133, (float).00180096606761054, (float)
		.00124753110589199, (float)9.22878876572938e-4, (float)
		7.14430421727287e-4, (float)5.71787281789705e-4, (float)
		4.69431007606482e-4, (float)3.93232835462917e-4, (float)
		3.34818889318298e-4, (float)2.88952148495752e-4, (float)
		2.52211615549573e-4, (float)2.22280580798883e-4, (float)
		1.97541838033063e-4, (float)1.76836855019718e-4, (float)
		1.59316899661821e-4, (float)1.44347930197334e-4, (float)
		1.31448068119965e-4, (float)1.20245444949303e-4, (float)
		1.10449144504599e-4, (float)1.01828770740567e-4, (float)
		9.41998224204238e-5, (float)8.74130545753834e-5, (float)
		8.13466262162801e-5, (float)7.59002269646219e-5, (float)
		-.00149282953213429, (float)-8.78204709546389e-4, (float)
		-5.02916549572035e-4, (float)-2.94822138512746e-4, (float)
		-1.75463996970783e-4, (float)-1.04008550460816e-4, (float)
		-5.96141953046458e-5, (float)-3.12038929076098e-5, (float)
		-1.2608973598023e-5, (float)-2.4289260857573e-7, (float)
		8.05996165414274e-6, (float)1.36507009262147e-5, (float)
		1.73964125472926e-5, (float)1.98672978842134e-5, (float)
		2.14463263790823e-5, (float)2.23954659232457e-5, (float)
		2.28967783814713e-5, (float)2.30785389811178e-5, (float)
		2.30321976080909e-5, (float)2.28236073720349e-5, (float)
		2.25005881105292e-5, (float)2.20981015361991e-5, (float)
		2.16418427448104e-5, (float)2.11507649256221e-5, (float)
		2.06388749782171e-5, (float)2.01165241997082e-5, (float)
		5.52213076721293e-4, (float)4.47932581552385e-4, (float)
		2.79520653992021e-4, (float)1.52468156198447e-4, (float)
		6.93271105657044e-5, (float)1.76258683069991e-5, (float)
		-1.35744996343269e-5, (float)-3.17972413350427e-5, (float)
		-4.18861861696693e-5, (float)-4.69004889379141e-5, (float)
		-4.87665447413787e-5, (float)-4.87010031186735e-5, (float)
		-4.74755620890087e-5, (float)-4.55813058138628e-5, (float)
		-4.33309644511266e-5, (float)-4.0923019315775e-5, (float)
		-3.84822638603221e-5, (float)-3.60857167535411e-5, (float)
		-3.37793306123367e-5, (float)-3.1588856077211e-5, (float)
		-2.95269561750807e-5, (float)-2.75978914828336e-5, (float)
		-2.58006174666884e-5, (float)-2.4130835676128e-5, (float)
		-2.25823509518346e-5, (float)-2.11479656768913e-5, (float)
		-4.7461779655996e-4, (float)-4.77864567147321e-4, (float)
		-3.20390228067038e-4, (float)-1.61105016119962e-4, (float)
		-4.25778101285435e-5, (float)3.44571294294968e-5, (float)
		7.97092684075675e-5, (float)1.03138236708272e-4, (float)
		1.12466775262204e-4, (float)1.13103642108481e-4, (float)
		1.08651634848774e-4, (float)1.01437951597662e-4, (float)
		9.29298396593364e-5, (float)8.4029313301609e-5, (float)
		7.52727991349134e-5, (float)6.69632521975731e-5, (float)
		5.92564547323195e-5, (float)5.22169308826976e-5, (float)
		4.58539485165361e-5, (float)4.01445513891487e-5, (float)
		3.50481730031328e-5, (float)3.05157995034347e-5, (float)
		2.64956119950516e-5, (float)2.29363633690998e-5, (float)
		1.97893056664022e-5, (float)1.70091984636413e-5, (float)
		7.36465810572578e-4, (float)8.72790805146194e-4, (float)
		6.22614862573135e-4, (float)2.85998154194304e-4, (float)
		3.84737672879366e-6, (float)-1.87906003636972e-4, (float)
		-2.97603646594555e-4, (float)-3.45998126832656e-4, (float)
		-3.53382470916038e-4, (float)-3.35715635775049e-4, (float)
		-3.0432112478904e-4, (float)-2.66722723047613e-4, (float)
		-2.2765421412282e-4, (float)-1.89922611854562e-4, (float)
		-1.55058918599094e-4, (float)-1.23778240761874e-4, (float)
		-9.62926147717644e-5, (float)-7.25178327714425e-5, (float)
		-5.22070028895634e-5, (float)-3.50347750511901e-5, (float)
		-2.06489761035552e-5, (float)-8.70106096849767e-6, (float)
		1.136986866751e-6, (float)9.16426474122779e-6, (float)
		1.56477785428873e-5, (float)2.08223629482467e-5 };


    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(), pow_dd(), log(), atan();

    /* Local variables */
#define alfa ((real *)&equiv_1)
#define beta ((real *)&equiv_4)
    static real relb, elim, rden;
    static integer kmax[5];
    static real crz32, asum, bsum, suma, sumb, upol[10];
#define alfa1 ((real *)&equiv_1)
#define alfa2 ((real *)&equiv_1 + 52)
#define beta1 ((real *)&equiv_4)
#define beta2 ((real *)&equiv_4 + 52)
#define beta3 ((real *)&equiv_4 + 104)
    static integer i__, j, k, l;
    static real z__;
    static integer iseta, isetb, klast;
    static real rzden;
    extern integer i1mach_();
    static real s1, t2;
    extern doublereal r1mach_();
    static integer kb;
    static real fi, ap, cr[10], dr[10];
    static integer jn;
    static real fn, sa, az;
    static integer jr;
    static real sb;
    static integer ks, ju, lr;
    static real ta, tb, z32, xx;
    static integer kstemp;
    static real fn2;
    static integer kp1;
    static real dfi, akm, phi, tfn, tau, rcz, tol, rtz, abw2, rfn2;
    static integer ksp1, lrp1;

/* ***BEGIN PROLOGUE  ASYJY */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BESJ and BESY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (ASYJY-S, DASYJY-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*                 ASYJY computes Bessel functions J and Y */
/*               for arguments X.GT.0.0 and orders FNU.GE.35.0 */
/*               on FLGJY = 1 and FLGJY = -1 respectively */

/*                                  INPUT */

/*      FUNJY - external function JAIRY or YAIRY */
/*          X - argument, X.GT.0.0E0 */
/*        FNU - order of the first Bessel function */
/*      FLGJY - selection flag */
/*              FLGJY =  1.0E0 gives the J function */
/*              FLGJY = -1.0E0 gives the Y function */
/*         IN - number of functions desired, IN = 1 or 2 */

/*                                  OUTPUT */

/*         Y  - a vector whose first in components contain the sequence */
/*       IFLW - a flag indicating underflow or overflow */
/*                    return variables for BESJ only */
/*      WK(1) = 1 - (X/FNU)**2 = W**2 */
/*      WK(2) = SQRT(ABS(WK(1))) */
/*      WK(3) = ABS(WK(2) - ATAN(WK(2)))  or */
/*              ABS(LN((1 + WK(2))/(X/FNU)) - WK(2)) */
/*            = ABS((2/3)*ZETA**(3/2)) */
/*      WK(4) = FNU*WK(3) */
/*      WK(5) = (1.5*WK(3)*FNU)**(1/3) = SQRT(ZETA)*FNU**(1/3) */
/*      WK(6) = SIGN(1.,W**2)*WK(5)**2 = SIGN(1.,W**2)*ZETA*FNU**(2/3) */
/*      WK(7) = FNU**(1/3) */

/*     Abstract */
/*         ASYJY implements the uniform asymptotic expansion of */
/*         the J and Y Bessel functions for FNU.GE.35 and real */
/*         X.GT.0.0E0. The forms are identical except for a change */
/*         in sign of some of the terms. This change in sign is */
/*         accomplished by means of the flag FLGJY = 1 or -1. On */
/*         FLGJY = 1 the AIRY functions AI(X) and DAI(X) are */
/*         supplied by the external function JAIRY, and on */
/*         FLGJY = -1 the AIRY functions BI(X) and DBI(X) are */
/*         supplied by the external function YAIRY. */

/* ***SEE ALSO  BESJ, BESY */
/* ***ROUTINES CALLED  I1MACH, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR section.  (WRB) */
/* ***END PROLOGUE  ASYJY */
    /* Parameter adjustments */
    --wk;
    --y;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  ASYJY */
    ta = r1mach_(&c__3);
    tol = dmax(ta,(float)1e-15);
    tb = r1mach_(&c__5);
    ju = i1mach_(&c__12);
    if (*flgjy == (float)1.) {
	goto L6;
    }
    jr = i1mach_(&c__11);
    elim = tb * (float)-2.303 * (ju + jr);
    goto L7;
L6:
    elim = (tb * ju + (float)3.) * (float)-2.303;
L7:
    fn = *fnu;
    *iflw = 0;
    i__1 = *in;
    for (jn = 1; jn <= i__1; ++jn) {
	xx = *x / fn;
	wk[1] = (float)1. - xx * xx;
	abw2 = dabs(wk[1]);
	wk[2] = sqrt(abw2);
	d__1 = (doublereal) fn;
	d__2 = (doublereal) con2;
	wk[7] = pow_dd(&d__1, &d__2);
	if (abw2 > (float).2775) {
	    goto L80;
	}

/*     ASYMPTOTIC EXPANSION */
/*     CASES NEAR X=FN, ABS(1.-(X/FN)**2).LE.0.2775 */
/*     COEFFICIENTS OF ASYMPTOTIC EXPANSION BY SERIES */

/*     ZETA AND TRUNCATION FOR A(ZETA) AND B(ZETA) SERIES */

/*     KMAX IS TRUNCATION INDEX FOR A(ZETA) AND B(ZETA) SERIES=MAX(2,S
A) */

	sa = (float)0.;
	if (abw2 == (float)0.) {
	    goto L10;
	}
	sa = tols / log(abw2);
L10:
	sb = sa;
	for (i__ = 1; i__ <= 5; ++i__) {
	    akm = dmax(sa,(float)2.);
	    kmax[i__ - 1] = (integer) akm;
	    sa += sb;
/* L20: */
	}
	kb = kmax[4];
	klast = kb - 1;
	sa = gama[kb - 1];
	i__2 = klast;
	for (k = 1; k <= i__2; ++k) {
	    --kb;
	    sa = sa * wk[1] + gama[kb - 1];
/* L30: */
	}
	z__ = wk[1] * sa;
	az = dabs(z__);
	rtz = sqrt(az);
	wk[3] = con1 * az * rtz;
	wk[4] = wk[3] * fn;
	wk[5] = rtz * wk[7];
	wk[6] = -wk[5] * wk[5];
	if (z__ <= (float)0.) {
	    goto L35;
	}
	if (wk[4] > elim) {
	    goto L75;
	}
	wk[6] = -wk[6];
L35:
	phi = sqrt(sqrt(sa + sa + sa + sa));

/*     B(ZETA) FOR S=0 */

	kb = kmax[4];
	klast = kb - 1;
	sb = beta[kb - 1];
	i__2 = klast;
	for (k = 1; k <= i__2; ++k) {
	    --kb;
	    sb = sb * wk[1] + beta[kb - 1];
/* L40: */
	}
	ksp1 = 1;
	fn2 = fn * fn;
	rfn2 = (float)1. / fn2;
	rden = (float)1.;
	asum = (float)1.;
	relb = tol * dabs(sb);
	bsum = sb;
	for (ks = 1; ks <= 4; ++ks) {
	    ++ksp1;
	    rden *= rfn2;

/*     A(ZETA) AND B(ZETA) FOR S=1,2,3,4 */

	    kstemp = 5 - ks;
	    kb = kmax[kstemp - 1];
	    klast = kb - 1;
	    sa = alfa[kb + ks * 26 - 27];
	    sb = beta[kb + ksp1 * 26 - 27];
	    i__2 = klast;
	    for (k = 1; k <= i__2; ++k) {
		--kb;
		sa = sa * wk[1] + alfa[kb + ks * 26 - 27];
		sb = sb * wk[1] + beta[kb + ksp1 * 26 - 27];
/* L50: */
	    }
	    ta = sa * rden;
	    tb = sb * rden;
	    asum += ta;
	    bsum += tb;
	    if (dabs(ta) <= tol && dabs(tb) <= relb) {
		goto L70;
	    }
/* L60: */
	}
L70:
	bsum /= fn * wk[7];
	goto L160;

L75:
	*iflw = 1;
	return 0;

L80:
	upol[0] = (float)1.;
	tau = (float)1. / wk[2];
	t2 = (float)1. / wk[1];
	if (wk[1] >= (float)0.) {
	    goto L90;
	}

/*     CASES FOR (X/FN).GT.SQRT(1.2775) */

	wk[3] = (r__1 = wk[2] - atan(wk[2]), dabs(r__1));
	wk[4] = wk[3] * fn;
	rcz = -con1 / wk[4];
	z32 = wk[3] * (float)1.5;
	d__1 = (doublereal) z32;
	d__2 = (doublereal) con2;
	rtz = pow_dd(&d__1, &d__2);
	wk[5] = rtz * wk[7];
	wk[6] = -wk[5] * wk[5];
	goto L100;
L90:

/*     CASES FOR (X/FN).LT.SQRT(0.7225) */

	wk[3] = (r__1 = log((wk[2] + (float)1.) / xx) - wk[2], dabs(r__1));
	wk[4] = wk[3] * fn;
	rcz = con1 / wk[4];
	if (wk[4] > elim) {
	    goto L75;
	}
	z32 = wk[3] * (float)1.5;
	d__1 = (doublereal) z32;
	d__2 = (doublereal) con2;
	rtz = pow_dd(&d__1, &d__2);
	d__1 = (doublereal) fn;
	d__2 = (doublereal) con2;
	wk[7] = pow_dd(&d__1, &d__2);
	wk[5] = rtz * wk[7];
	wk[6] = wk[5] * wk[5];
L100:
	phi = sqrt((rtz + rtz) * tau);
	tb = (float)1.;
	asum = (float)1.;
	tfn = tau / fn;
	rden = (float)1. / fn;
	rfn2 = rden * rden;
	rden = (float)1.;
	upol[1] = (c__[0] * t2 + c__[1]) * tfn;
	crz32 = con548 * rcz;
	bsum = upol[1] + crz32;
	relb = tol * dabs(bsum);
	ap = tfn;
	ks = 0;
	kp1 = 2;
	rzden = rcz;
	l = 2;
	iseta = 0;
	isetb = 0;
	for (lr = 2; lr <= 8; lr += 2) {

/*     COMPUTE TWO U POLYNOMIALS FOR NEXT A(ZETA) AND B(ZETA) */

	    lrp1 = lr + 1;
	    i__2 = lrp1;
	    for (k = lr; k <= i__2; ++k) {
		++ks;
		++kp1;
		++l;
		s1 = c__[l - 1];
		i__3 = kp1;
		for (j = 2; j <= i__3; ++j) {
		    ++l;
		    s1 = s1 * t2 + c__[l - 1];
/* L110: */
		}
		ap *= tfn;
		upol[kp1 - 1] = ap * s1;
		cr[ks - 1] = br[ks - 1] * rzden;
		rzden *= rcz;
		dr[ks - 1] = ar[ks - 1] * rzden;
/* L120: */
	    }
	    suma = upol[lrp1 - 1];
	    sumb = upol[lr + 1] + upol[lrp1 - 1] * crz32;
	    ju = lrp1;
	    i__2 = lr;
	    for (jr = 1; jr <= i__2; ++jr) {
		--ju;
		suma += cr[jr - 1] * upol[ju - 1];
		sumb += dr[jr - 1] * upol[ju - 1];
/* L130: */
	    }
	    rden *= rfn2;
	    tb = -tb;
	    if (wk[1] > (float)0.) {
		tb = dabs(tb);
	    }
	    if (rden < tol) {
		goto L131;
	    }
	    asum += suma * tb;
	    bsum += sumb * tb;
	    goto L140;
L131:
	    if (iseta == 1) {
		goto L132;
	    }
	    if (dabs(suma) < tol) {
		iseta = 1;
	    }
	    asum += suma * tb;
L132:
	    if (isetb == 1) {
		goto L133;
	    }
	    if (dabs(sumb) < relb) {
		isetb = 1;
	    }
	    bsum += sumb * tb;
L133:
	    if (iseta == 1 && isetb == 1) {
		goto L150;
	    }
L140:
	    ;
	}
L150:
	tb = wk[5];
	if (wk[1] > (float)0.) {
	    tb = -tb;
	}
	bsum /= tb;

L160:
	(*funjy)(&wk[6], &wk[5], &wk[4], &fi, &dfi);
	ta = (float)1. / tol;
	tb = r1mach_(&c__1) * ta * (float)1e3;
	if (dabs(fi) > tb) {
	    goto L165;
	}
	fi *= ta;
	dfi *= ta;
	phi *= tol;
L165:
	y[jn] = *flgjy * phi * (fi * asum + dfi * bsum) / wk[7];
	fn -= *flgjy;
/* L170: */
    }
    return 0;
} /* asyjy_ */

#undef beta3
#undef beta2
#undef beta1
#undef alfa2
#undef alfa1
#undef beta
#undef alfa


