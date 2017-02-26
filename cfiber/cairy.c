/* cairy.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__5 = 5;
static integer c__11 = 11;
static integer c__9 = 9;
static integer c__1 = 1;

/* DECK CAIRY */
/* Subroutine */ int cairy_(z__, id, kode, ai, nz, ierr)
complex *z__;
integer *id, *kode;
complex *ai;
integer *nz, *ierr;
{
    /* Initialized data */

    static real tth = (float).666666666666666667;
    static real c1 = (float).35502805388781724;
    static real c2 = (float).258819403792806799;
    static real coef = (float).183776298473930683;
    static complex cone = {(float)1.,(float)0.};

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    doublereal d__1, d__2;
    complex q__1, q__2, q__3, q__4, q__5, q__6;

    /* Builtin functions */
    double c_abs(), r_imag();
    void c_sqrt(), c_exp();
    double log(), pow_dd(), sqrt();

    /* Local variables */
    static real sfac, alim, elim, alaz, atrm;
    extern /* Subroutine */ int cacai_();
    static integer k, iflag;
    extern /* Subroutine */ int cbknu_();
    static real d1, d2;
    static integer k1, k2;
    extern integer i1mach_();
    static complex s1, s2;
    extern doublereal r1mach_();
    static complex z3;
    static real aa, bb, ad, ak, bk, ck, dk, az;
    static complex cy[1];
    static integer nn;
    static real rl;
    static integer mr;
    static real zi, zr, az3, z3i, z3r, fid, dig, r1m5;
    static complex csq;
    static real fnu;
    static complex zta;
    static real tol;
    static complex trm1, trm2;

/* ***BEGIN PROLOGUE  CAIRY */
/* ***PURPOSE  Compute the Airy function Ai(z) or its derivative dAi/dz */
/*            for complex argument z.  A scaling option is available */
/*            to help avoid underflow and overflow. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10D */
/* ***TYPE      COMPLEX (CAIRY-C, ZAIRY-C) */
/* ***KEYWORDS  AIRY FUNCTION, BESSEL FUNCTION OF ORDER ONE THIRD, */
/*             BESSEL FUNCTION OF ORDER TWO THIRDS */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*         On KODE=1, CAIRY computes the complex Airy function Ai(z) */
/*         or its derivative dAi/dz on ID=0 or ID=1 respectively. On */
/*         KODE=2, a scaling option exp(zeta)*Ai(z) or exp(zeta)*dAi/dz */
/*         is provided to remove the exponential decay in -pi/3<arg(z) */
/*         <pi/3 and the exponential growth in pi/3<abs(arg(z))<pi where 
*/
/*         zeta=(2/3)*z**(3/2). */

/*         While the Airy functions Ai(z) and dAi/dz are analytic in */
/*         the whole z-plane, the corresponding scaled functions defined 
*/
/*         for KODE=2 have a cut along the negative real axis. */

/*         Input */
/*           Z      - Argument of type COMPLEX */
/*           ID     - Order of derivative, ID=0 or ID=1 */
/*           KODE   - A parameter to indicate the scaling option */
/*                    KODE=1  returns */
/*                            AI=Ai(z)  on ID=0 */
/*                            AI=dAi/dz on ID=1 */
/*                            at z=Z */
/*                        =2  returns */
/*                            AI=exp(zeta)*Ai(z)  on ID=0 */
/*                            AI=exp(zeta)*dAi/dz on ID=1 */
/*                            at z=Z where zeta=(2/3)*z**(3/2) */

/*         Output */
/*           AI     - Result of type COMPLEX */
/*           NZ     - Underflow indicator */
/*                    NZ=0    Normal return */
/*                    NZ=1    AI=0 due to underflow in */
/*                            -pi/3<arg(Z)<pi/3 on KODE=1 */
/*           IERR   - Error flag */
/*                    IERR=0  Normal return     - COMPUTATION COMPLETED */
/*                    IERR=1  Input error       - NO COMPUTATION */
/*                    IERR=2  Overflow          - NO COMPUTATION */
/*                            (Re(Z) too large with KODE=1) */
/*                    IERR=3  Precision warning - COMPUTATION COMPLETED */
/*                            (Result has less than half precision) */
/*                    IERR=4  Precision error   - NO COMPUTATION */
/*                            (Result has no precision) */
/*                    IERR=5  Algorithmic error - NO COMPUTATION */
/*                            (Termination condition not met) */

/* *Long Description: */

/*         Ai(z) and dAi/dz are computed from K Bessel functions by */

/*                Ai(z) =  c*sqrt(z)*K(1/3,zeta) */
/*               dAi/dz = -c*   z   *K(2/3,zeta) */
/*                    c =  1/(pi*sqrt(3)) */
/*                 zeta =  (2/3)*z**(3/2) */

/*         when abs(z)>1 and from power series when abs(z)<=1. */

/*         In most complex variable computation, one must evaluate ele- */
/*         mentary functions.  When the magnitude of Z is large, losses */
/*         of significance by argument reduction occur.  Consequently, if 
*/
/*         the magnitude of ZETA=(2/3)*Z**(3/2) exceeds U1=SQRT(0.5/UR), 
*/
/*         then losses exceeding half precision are likely and an error */
/*         flag IERR=3 is triggered where UR=R1MACH(4)=UNIT ROUNDOFF. */
/*         Also, if the magnitude of ZETA is larger than U2=0.5/UR, then 
*/
/*         all significance is lost and IERR=4.  In order to use the INT 
*/
/*         function, ZETA must be further restricted not to exceed */
/*         U3=I1MACH(9)=LARGEST INTEGER.  Thus, the magnitude of ZETA */
/*         must be restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, 
*/
/*         and U3 are approximately 2.0E+3, 4.2E+6, 2.1E+9 in single */
/*         precision and 4.7E+7, 2.3E+15, 2.1E+9 in double precision. */
/*         This makes U2 limiting is single precision and U3 limiting */
/*         in double precision.  This means that the magnitude of Z */
/*         cannot exceed approximately 3.4E+4 in single precision and */
/*         2.1E+6 in double precision.  This also means that one can */
/*         expect to retain, in the worst cases on 32-bit machines, */
/*         no digits in single precision and only 6 digits in double */
/*         precision. */

/*         The approximate relative error in the magnitude of a complex */
/*         Bessel function can be expressed as P*10**S where P=MAX(UNIT */
/*         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre- */
/*         sents the increase in error due to argument reduction in the */
/*         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))), */
/*         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF 
*/
/*         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may 
*/
/*         have only absolute accuracy.  This is most likely to occur */
/*         when one component (in magnitude) is larger than the other by 
*/
/*         several orders of magnitude.  If one component is 10**K larger 
*/
/*         than the other, then one can expect only MAX(ABS(LOG10(P))-K, 
*/
/*         0) significant digits; or, stated another way, when K exceeds 
*/
/*         the exponent of P, no significant digits remain in the smaller 
*/
/*         component.  However, the phase angle retains absolute accuracy 
*/
/*         because, in complex arithmetic with precision P, the smaller */
/*         component will not (as a rule) decrease below P times the */
/*         magnitude of the larger component. In these extreme cases, */
/*         the principal phase angle is on the order of +P, -P, PI/2-P, */
/*         or -PI/2+P. */

/* ***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe- */
/*                 matical Functions, National Bureau of Standards */
/*                 Applied Mathematics Series 55, U. S. Department */
/*                 of Commerce, Tenth Printing (1972) or later. */
/*               2. D. E. Amos, Computation of Bessel Functions of */
/*                 Complex Argument and Large Order, Report SAND83-0643, 
*/
/*                 Sandia National Laboratories, Albuquerque, NM, May */
/*                 1983. */
/*               3. D. E. Amos, A Subroutine Package for Bessel Functions 
*/
/*                 of a Complex Argument and Nonnegative Order, Report */
/*                 SAND85-1018, Sandia National Laboratory, Albuquerque, 
*/
/*                 NM, May 1985. */
/*               4. D. E. Amos, A portable package for Bessel functions */
/*                 of a complex argument and nonnegative order, ACM */
/*                 Transactions on Mathematical Software, 12 (September */
/*                 1986), pp. 265-273. */

/* ***ROUTINES CALLED  CACAI, CBKNU, I1MACH, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   890801  REVISION DATE from Version 3.2 */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   920128  Category corrected.  (WRB) */
/*   920811  Prologue revised.  (DWL) */
/* ***END PROLOGUE  CAIRY */
/* ***FIRST EXECUTABLE STATEMENT  CAIRY */
    *ierr = 0;
    *nz = 0;
    if (*id < 0 || *id > 1) {
	*ierr = 1;
    }
    if (*kode < 1 || *kode > 2) {
	*ierr = 1;
    }
    if (*ierr != 0) {
	return 0;
    }
    az = c_abs(z__);
/* Computing MAX */
    r__1 = r1mach_(&c__4);
    tol = dmax(r__1,(float)1e-18);
    fid = (real) (*id);
    if (az > (float)1.) {
	goto L60;
    }
/* -----------------------------------------------------------------------
 */
/*     POWER SERIES FOR ABS(Z).LE.1. */
/* -----------------------------------------------------------------------
 */
    s1.r = cone.r, s1.i = cone.i;
    s2.r = cone.r, s2.i = cone.i;
    if (az < tol) {
	goto L160;
    }
    aa = az * az;
    if (aa < tol / az) {
	goto L40;
    }
    trm1.r = cone.r, trm1.i = cone.i;
    trm2.r = cone.r, trm2.i = cone.i;
    atrm = (float)1.;
    q__2.r = z__->r * z__->r - z__->i * z__->i, q__2.i = z__->r * z__->i + 
	    z__->i * z__->r;
    q__1.r = q__2.r * z__->r - q__2.i * z__->i, q__1.i = q__2.r * z__->i + 
	    q__2.i * z__->r;
    z3.r = q__1.r, z3.i = q__1.i;
    az3 = az * aa;
    ak = fid + (float)2.;
    bk = (float)3. - fid - fid;
    ck = (float)4. - fid;
    dk = fid + (float)3. + fid;
    d1 = ak * dk;
    d2 = bk * ck;
    ad = dmin(d1,d2);
    ak = fid * (float)9. + (float)24.;
    bk = (float)30. - fid * (float)9.;
    z3r = z3.r;
    z3i = r_imag(&z3);
    for (k = 1; k <= 25; ++k) {
	r__1 = z3r / d1;
	r__2 = z3i / d1;
	q__2.r = r__1, q__2.i = r__2;
	q__1.r = trm1.r * q__2.r - trm1.i * q__2.i, q__1.i = trm1.r * q__2.i 
		+ trm1.i * q__2.r;
	trm1.r = q__1.r, trm1.i = q__1.i;
	q__1.r = s1.r + trm1.r, q__1.i = s1.i + trm1.i;
	s1.r = q__1.r, s1.i = q__1.i;
	r__1 = z3r / d2;
	r__2 = z3i / d2;
	q__2.r = r__1, q__2.i = r__2;
	q__1.r = trm2.r * q__2.r - trm2.i * q__2.i, q__1.i = trm2.r * q__2.i 
		+ trm2.i * q__2.r;
	trm2.r = q__1.r, trm2.i = q__1.i;
	q__1.r = s2.r + trm2.r, q__1.i = s2.i + trm2.i;
	s2.r = q__1.r, s2.i = q__1.i;
	atrm = atrm * az3 / ad;
	d1 += ak;
	d2 += bk;
	ad = dmin(d1,d2);
	if (atrm < tol * ad) {
	    goto L40;
	}
	ak += (float)18.;
	bk += (float)18.;
/* L30: */
    }
L40:
    if (*id == 1) {
	goto L50;
    }
    q__3.r = c1, q__3.i = (float)0.;
    q__2.r = s1.r * q__3.r - s1.i * q__3.i, q__2.i = s1.r * q__3.i + s1.i * 
	    q__3.r;
    q__5.r = z__->r * s2.r - z__->i * s2.i, q__5.i = z__->r * s2.i + z__->i * 
	    s2.r;
    q__6.r = c2, q__6.i = (float)0.;
    q__4.r = q__5.r * q__6.r - q__5.i * q__6.i, q__4.i = q__5.r * q__6.i + 
	    q__5.i * q__6.r;
    q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
    ai->r = q__1.r, ai->i = q__1.i;
    if (*kode == 1) {
	return 0;
    }
    c_sqrt(&q__3, z__);
    q__2.r = z__->r * q__3.r - z__->i * q__3.i, q__2.i = z__->r * q__3.i + 
	    z__->i * q__3.r;
    q__4.r = tth, q__4.i = (float)0.;
    q__1.r = q__2.r * q__4.r - q__2.i * q__4.i, q__1.i = q__2.r * q__4.i + 
	    q__2.i * q__4.r;
    zta.r = q__1.r, zta.i = q__1.i;
    c_exp(&q__2, &zta);
    q__1.r = ai->r * q__2.r - ai->i * q__2.i, q__1.i = ai->r * q__2.i + ai->i 
	    * q__2.r;
    ai->r = q__1.r, ai->i = q__1.i;
    return 0;
L50:
    q__2.r = -s2.r, q__2.i = -s2.i;
    q__3.r = c2, q__3.i = (float)0.;
    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * q__3.i + 
	    q__2.i * q__3.r;
    ai->r = q__1.r, ai->i = q__1.i;
    if (az > tol) {
	q__4.r = z__->r * z__->r - z__->i * z__->i, q__4.i = z__->r * z__->i 
		+ z__->i * z__->r;
	q__3.r = q__4.r * s1.r - q__4.i * s1.i, q__3.i = q__4.r * s1.i + 
		q__4.i * s1.r;
	r__1 = c1 / (fid + (float)1.);
	q__5.r = r__1, q__5.i = (float)0.;
	q__2.r = q__3.r * q__5.r - q__3.i * q__5.i, q__2.i = q__3.r * q__5.i 
		+ q__3.i * q__5.r;
	q__1.r = ai->r + q__2.r, q__1.i = ai->i + q__2.i;
	ai->r = q__1.r, ai->i = q__1.i;
    }
    if (*kode == 1) {
	return 0;
    }
    c_sqrt(&q__3, z__);
    q__2.r = z__->r * q__3.r - z__->i * q__3.i, q__2.i = z__->r * q__3.i + 
	    z__->i * q__3.r;
    q__4.r = tth, q__4.i = (float)0.;
    q__1.r = q__2.r * q__4.r - q__2.i * q__4.i, q__1.i = q__2.r * q__4.i + 
	    q__2.i * q__4.r;
    zta.r = q__1.r, zta.i = q__1.i;
    c_exp(&q__2, &zta);
    q__1.r = ai->r * q__2.r - ai->i * q__2.i, q__1.i = ai->r * q__2.i + ai->i 
	    * q__2.r;
    ai->r = q__1.r, ai->i = q__1.i;
    return 0;
/* -----------------------------------------------------------------------
 */
/*     CASE FOR ABS(Z).GT.1.0 */
/* -----------------------------------------------------------------------
 */
L60:
    fnu = (fid + (float)1.) / (float)3.;
/* -----------------------------------------------------------------------
 */
/*     SET PARAMETERS RELATED TO MACHINE CONSTANTS. */
/*     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18. */
/*     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. */
/*     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND */
/*     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR */
/*     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. */
/*     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. 
*/
/*     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). */
/* -----------------------------------------------------------------------
 */
    k1 = i1mach_(&c__12);
    k2 = i1mach_(&c__13);
    r1m5 = r1mach_(&c__5);
/* Computing MIN */
    i__1 = abs(k1), i__2 = abs(k2);
    k = min(i__1,i__2);
    elim = (k * r1m5 - (float)3.) * (float)2.303;
    k1 = i1mach_(&c__11) - 1;
    aa = r1m5 * k1;
    dig = dmin(aa,(float)18.);
    aa *= (float)2.303;
/* Computing MAX */
    r__1 = -aa;
    alim = elim + dmax(r__1,(float)-41.45);
    rl = dig * (float)1.2 + (float)3.;
    alaz = log(az);
/* -----------------------------------------------------------------------
 */
/*     TEST FOR RANGE */
/* -----------------------------------------------------------------------
 */
    aa = (float).5 / tol;
    bb = i1mach_(&c__9) * (float).5;
    aa = dmin(aa,bb);
    d__1 = (doublereal) aa;
    d__2 = (doublereal) tth;
    aa = pow_dd(&d__1, &d__2);
    if (az > aa) {
	goto L260;
    }
    aa = sqrt(aa);
    if (az > aa) {
	*ierr = 3;
    }
    c_sqrt(&q__1, z__);
    csq.r = q__1.r, csq.i = q__1.i;
    q__2.r = z__->r * csq.r - z__->i * csq.i, q__2.i = z__->r * csq.i + 
	    z__->i * csq.r;
    q__3.r = tth, q__3.i = (float)0.;
    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * q__3.i + 
	    q__2.i * q__3.r;
    zta.r = q__1.r, zta.i = q__1.i;
/* -----------------------------------------------------------------------
 */
/*     RE(ZTA).LE.0 WHEN RE(Z).LT.0, ESPECIALLY WHEN IM(Z) IS SMALL */
/* -----------------------------------------------------------------------
 */
    iflag = 0;
    sfac = (float)1.;
    zi = r_imag(z__);
    zr = z__->r;
    ak = r_imag(&zta);
    if (zr >= (float)0.) {
	goto L70;
    }
    bk = zta.r;
    ck = -dabs(bk);
    q__1.r = ck, q__1.i = ak;
    zta.r = q__1.r, zta.i = q__1.i;
L70:
    if (zi != (float)0.) {
	goto L80;
    }
    if (zr > (float)0.) {
	goto L80;
    }
    q__1.r = (float)0., q__1.i = ak;
    zta.r = q__1.r, zta.i = q__1.i;
L80:
    aa = zta.r;
    if (aa >= (float)0. && zr > (float)0.) {
	goto L100;
    }
    if (*kode == 2) {
	goto L90;
    }
/* -----------------------------------------------------------------------
 */
/*     OVERFLOW TEST */
/* -----------------------------------------------------------------------
 */
    if (aa > -alim) {
	goto L90;
    }
    aa = -aa + alaz * (float).25;
    iflag = 1;
    sfac = tol;
    if (aa > elim) {
	goto L240;
    }
L90:
/* -----------------------------------------------------------------------
 */
/*     CBKNU AND CACAI RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2 */
/* -----------------------------------------------------------------------
 */
    mr = 1;
    if (zi < (float)0.) {
	mr = -1;
    }
    cacai_(&zta, &fnu, kode, &mr, &c__1, cy, &nn, &rl, &tol, &elim, &alim);
    if (nn < 0) {
	goto L250;
    }
    *nz += nn;
    goto L120;
L100:
    if (*kode == 2) {
	goto L110;
    }
/* -----------------------------------------------------------------------
 */
/*     UNDERFLOW TEST */
/* -----------------------------------------------------------------------
 */
    if (aa < alim) {
	goto L110;
    }
    aa = -aa - alaz * (float).25;
    iflag = 2;
    sfac = (float)1. / tol;
    if (aa < -elim) {
	goto L180;
    }
L110:
    cbknu_(&zta, &fnu, kode, &c__1, cy, nz, &tol, &elim, &alim);
L120:
    q__2.r = coef, q__2.i = (float)0.;
    q__1.r = cy[0].r * q__2.r - cy[0].i * q__2.i, q__1.i = cy[0].r * q__2.i + 
	    cy[0].i * q__2.r;
    s1.r = q__1.r, s1.i = q__1.i;
    if (iflag != 0) {
	goto L140;
    }
    if (*id == 1) {
	goto L130;
    }
    q__1.r = csq.r * s1.r - csq.i * s1.i, q__1.i = csq.r * s1.i + csq.i * 
	    s1.r;
    ai->r = q__1.r, ai->i = q__1.i;
    return 0;
L130:
    q__2.r = -z__->r, q__2.i = -z__->i;
    q__1.r = q__2.r * s1.r - q__2.i * s1.i, q__1.i = q__2.r * s1.i + q__2.i * 
	    s1.r;
    ai->r = q__1.r, ai->i = q__1.i;
    return 0;
L140:
    q__2.r = sfac, q__2.i = (float)0.;
    q__1.r = s1.r * q__2.r - s1.i * q__2.i, q__1.i = s1.r * q__2.i + s1.i * 
	    q__2.r;
    s1.r = q__1.r, s1.i = q__1.i;
    if (*id == 1) {
	goto L150;
    }
    q__1.r = s1.r * csq.r - s1.i * csq.i, q__1.i = s1.r * csq.i + s1.i * 
	    csq.r;
    s1.r = q__1.r, s1.i = q__1.i;
    r__1 = (float)1. / sfac;
    q__2.r = r__1, q__2.i = (float)0.;
    q__1.r = s1.r * q__2.r - s1.i * q__2.i, q__1.i = s1.r * q__2.i + s1.i * 
	    q__2.r;
    ai->r = q__1.r, ai->i = q__1.i;
    return 0;
L150:
    q__2.r = -s1.r, q__2.i = -s1.i;
    q__1.r = q__2.r * z__->r - q__2.i * z__->i, q__1.i = q__2.r * z__->i + 
	    q__2.i * z__->r;
    s1.r = q__1.r, s1.i = q__1.i;
    r__1 = (float)1. / sfac;
    q__2.r = r__1, q__2.i = (float)0.;
    q__1.r = s1.r * q__2.r - s1.i * q__2.i, q__1.i = s1.r * q__2.i + s1.i * 
	    q__2.r;
    ai->r = q__1.r, ai->i = q__1.i;
    return 0;
L160:
    aa = r1mach_(&c__1) * (float)1e3;
    s1.r = (float)0., s1.i = (float)0.;
    if (*id == 1) {
	goto L170;
    }
    if (az > aa) {
	q__2.r = c2, q__2.i = (float)0.;
	q__1.r = q__2.r * z__->r - q__2.i * z__->i, q__1.i = q__2.r * z__->i 
		+ q__2.i * z__->r;
	s1.r = q__1.r, s1.i = q__1.i;
    }
    q__2.r = c1, q__2.i = (float)0.;
    q__1.r = q__2.r - s1.r, q__1.i = q__2.i - s1.i;
    ai->r = q__1.r, ai->i = q__1.i;
    return 0;
L170:
    q__2.r = c2, q__2.i = (float)0.;
    q__1.r = -q__2.r, q__1.i = -q__2.i;
    ai->r = q__1.r, ai->i = q__1.i;
    aa = sqrt(aa);
    if (az > aa) {
	q__2.r = z__->r * z__->r - z__->i * z__->i, q__2.i = z__->r * z__->i 
		+ z__->i * z__->r;
	q__1.r = q__2.r * (float).5 - q__2.i * (float)0., q__1.i = q__2.r * (
		float)0. + q__2.i * (float).5;
	s1.r = q__1.r, s1.i = q__1.i;
    }
    q__3.r = c1, q__3.i = (float)0.;
    q__2.r = s1.r * q__3.r - s1.i * q__3.i, q__2.i = s1.r * q__3.i + s1.i * 
	    q__3.r;
    q__1.r = ai->r + q__2.r, q__1.i = ai->i + q__2.i;
    ai->r = q__1.r, ai->i = q__1.i;
    return 0;
L180:
    *nz = 1;
    ai->r = (float)0., ai->i = (float)0.;
    return 0;
L240:
    *nz = 0;
    *ierr = 2;
    return 0;
L250:
    if (nn == -1) {
	goto L240;
    }
    *nz = 0;
    *ierr = 5;
    return 0;
L260:
    *ierr = 4;
    *nz = 0;
    return 0;
} /* cairy_ */

