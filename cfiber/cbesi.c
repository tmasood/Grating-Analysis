/* cbesi.f -- translated by f2c (version 20000121).
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

/* DECK CBESI */
/* Subroutine */ int cbesi(complex *z__, real *fnu, integer *kode, integer *n,
			    complex *cy, integer *nz, integer *ierr)
{
    /* Initialized data */

    static real pi = (float)3.14159265358979324;
    static complex cone = {(float)1.,(float)0.};

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    complex q__1, q__2;

    /* Builtin functions */
    double r_imag(), c_abs(), sqrt(), cos(), sin();

    /* Local variables */
    static real alim, elim;
    static complex csgn;
    static real atol, fnul, rtol;
    static integer i__, k;
    static real ascle;
    extern /* Subroutine */ int cbinu_(complex *, real *, integer *,
				       integer *, complex *, integer *,
				       real *, real *, real *, real *, real *);
    static integer k1, k2;
    extern integer i1mach_();
    static real s1, s2;
    extern doublereal r1mach_();
    static real aa, bb, fn, az;
    static integer nn;
    static real rl;
    static complex zn;
    static real xx, yy, dig, arg, r1m5;
    static integer inu;
    static real tol;

/* ***BEGIN PROLOGUE  CBESI */
/* ***PURPOSE  Compute a sequence of the Bessel functions I(a,z) for */
/*            complex argument z and real nonnegative orders a=b,b+1, */
/*            b+2,... where b>0.  A scaling option is available to */
/*            help avoid overflow. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10B4 */
/* ***TYPE      COMPLEX (CBESI-C, ZBESI-C) */
/* ***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT, I BESSEL FUNCTIONS, */
/*             MODIFIED BESSEL FUNCTIONS */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*         On KODE=1, CBESI computes an N-member sequence of complex */
/*         Bessel functions CY(L)=I(FNU+L-1,Z) for real nonnegative */
/*         orders FNU+L-1, L=1,...,N and complex Z in the cut plane */
/*         -pi<arg(Z)<=pi.  On KODE=2, CBESI returns the scaled functions */

/*            CY(L) = exp(-abs(X))*I(FNU+L-1,Z), L=1,...,N and X=Re(Z) */

/*         which removes the exponential growth in both the left and */
/*         right half-planes as Z goes to infinity. */

/*         Input */
/*           Z      - Argument of type COMPLEX */
/*           FNU    - Initial order of type REAL, FNU>=0 */
/*           KODE   - A parameter to indicate the scaling option */
/*                    KODE=1  returns */
/*                            CY(L)=I(FNU+L-1,Z), L=1,...,N */
/*                        =2  returns */
/*                            CY(L)=exp(-abs(X))*I(FNU+L-1,Z), L=1,...,N */
/*                            where X=Re(Z) */
/*           N      - Number of terms in the sequence, N>=1 */

/*         Output */
/*           CY     - Result vector of type COMPLEX */
/*           NZ     - Number of underflows set to zero */
/*                    NZ=0    Normal return */
/*                    NZ>0    CY(L)=0, L=N-NZ+1,...,N */
/*           IERR   - Error flag */
/*                    IERR=0  Normal return     - COMPUTATION COMPLETED */
/*                    IERR=1  Input error       - NO COMPUTATION */
/*                    IERR=2  Overflow          - NO COMPUTATION */
/*                            (Re(Z) too large on KODE=1) */
/*                    IERR=3  Precision warning - COMPUTATION COMPLETED */
/*                            (Result has half precision or less */
/*                            because abs(Z) or FNU+N-1 is large) */
/*                    IERR=4  Precision error   - NO COMPUTATION */
/*                            (Result has no precision because */
/*                            abs(Z) or FNU+N-1 is too large) */
/*                    IERR=5  Algorithmic error - NO COMPUTATION */
/*                            (Termination condition not met) */

/* *Long Description: */

/*         The computation of I(a,z) is carried out by the power series */
/*         for small abs(z), the asymptotic expansion for large abs(z), */
/*         the Miller algorithm normalized by the Wronskian and a */
/*         Neumann series for intermediate magnitudes of z, and the */
/*         uniform asymptotic expansions for I(a,z) and J(a,z) for */
/*         large orders a.  Backward recurrence is used to generate */
/*         sequences or reduce orders when necessary. */

/*         The calculations above are done in the right half plane and */
/*         continued into the left half plane by the formula */

/*            I(a,z*exp(t)) = exp(t*a)*I(a,z), Re(z)>0 */
/*                        t = i*pi or -i*pi */

/*         For negative orders, the formula */

/*            I(-a,z) = I(a,z) + (2/pi)*sin(pi*a)*K(a,z) */

/*         can be used.  However, for large orders close to integers the */
/*         the function changes radically.  When a is a large positive */
/*         integer, the magnitude of I(-a,z)=I(a,z) is a large */
/*         negative power of ten. But when a is not an integer, */
/*         K(a,z) dominates in magnitude with a large positive power of */
/*         ten and the most that the second term can be reduced is by */
/*         unit roundoff from the coefficient. Thus, wide changes can */
/*         occur within unit roundoff of a large integer for a. Here, */
/*         large means a>abs(z). */

/*         In most complex variable computation, one must evaluate ele- */
/*         mentary functions.  When the magnitude of Z or FNU+N-1 is */
/*         large, losses of significance by argument reduction occur. */
/*         Consequently, if either one exceeds U1=SQRT(0.5/UR), then */
/*         losses exceeding half precision are likely and an error flag */
/*         IERR=3 is triggered where UR=R1MACH(4)=UNIT ROUNDOFF.  Also, */
/*         if either is larger than U2=0.5/UR, then all significance is */
/*         lost and IERR=4.  In order to use the INT function, arguments */
/*         must be further restricted not to exceed the largest machine */
/*         integer, U3=I1MACH(9).  Thus, the magnitude of Z and FNU+N-1 */
/*         is restricted by MIN(U2,U3).  In IEEE arithmetic, U1,U2, and */
/*         U3 approximate 2.0E+3, 4.2E+6, 2.1E+9 in single precision */
/*         and 4.7E+7, 2.3E+15 and 2.1E+9 in double precision.  This */
/*         makes U2 limiting in single precision and U3 limiting in */
/*         double precision.  This means that one can expect to retain, */
/*         in the worst cases on IEEE machines, no digits in single pre- */
/*         cision and only 6 digits in double precision.  Similar con- */
/*         siderations hold for other machines. */

/*         The approximate relative error in the magnitude of a complex */
/*         Bessel function can be expressed as P*10**S where P=MAX(UNIT */
/*         ROUNDOFF,1.0E-18) is the nominal precision and 10**S repre- */
/*         sents the increase in error due to argument reduction in the */
/*         elementary functions.  Here, S=MAX(1,ABS(LOG10(ABS(Z))), */
/*         ABS(LOG10(FNU))) approximately (i.e., S=MAX(1,ABS(EXPONENT OF */
/*         ABS(Z),ABS(EXPONENT OF FNU)) ).  However, the phase angle may */
/*         have only absolute accuracy.  This is most likely to occur */
/*         when one component (in magnitude) is larger than the other by */
/*         several orders of magnitude.  If one component is 10**K larger */
/*         than the other, then one can expect only MAX(ABS(LOG10(P))-K, */
/*         0) significant digits; or, stated another way, when K exceeds */
/*         the exponent of P, no significant digits remain in the smaller */
/*         component.  However, the phase angle retains absolute accuracy */
/*         because, in complex arithmetic with precision P, the smaller */
/*         component will not (as a rule) decrease below P times the */
/*         magnitude of the larger component.  In these extreme cases, */
/*         the principal phase angle is on the order of +P, -P, PI/2-P, */
/*         or -PI/2+P. */

/* ***REFERENCES  1. M. Abramowitz and I. A. Stegun, Handbook of Mathe- */
/*                 matical Functions, National Bureau of Standards */
/*                 Applied Mathematics Series 55, U. S. Department */
/*                 of Commerce, Tenth Printing (1972) or later. */
/*               2. D. E. Amos, Computation of Bessel Functions of */
/*                 Complex Argument, Report SAND83-0086, Sandia National */
/*                 Laboratories, Albuquerque, NM, May 1983. */
/*               3. D. E. Amos, Computation of Bessel Functions of */
/*                 Complex Argument and Large Order, Report SAND83-0643, */
/*                 Sandia National Laboratories, Albuquerque, NM, May */
/*                 1983. */
/*               4. D. E. Amos, A Subroutine Package for Bessel Functions */
/*                 of a Complex Argument and Nonnegative Order, Report */
/*                 SAND85-1018, Sandia National Laboratory, Albuquerque, */
/*                 NM, May 1985. */
/*               5. D. E. Amos, A portable package for Bessel functions */
/*                 of a complex argument and nonnegative order, ACM */
/*                 Transactions on Mathematical Software, 12 (September */
/*                 1986), pp. 265-273. */

/* ***ROUTINES CALLED  CBINU, I1MACH, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   890801  REVISION DATE from Version 3.2 */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   920128  Category corrected.  (WRB) */
/*   920811  Prologue revised.  (DWL) */
/* ***END PROLOGUE  CBESI */
    /* Parameter adjustments */
    --cy;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  CBESI */
    *ierr = 0;
    *nz = 0;
    if (*fnu < (float)0.) {
	*ierr = 1;
    }
    if (*kode < 1 || *kode > 2) {
	*ierr = 1;
    }
    if (*n < 1) {
	*ierr = 1;
    }
    if (*ierr != 0) {
	return 0;
    }
    xx = z__->r;
    yy = r_imag(z__);
/* ----------------------------------------------------------------------- */
/*     SET PARAMETERS RELATED TO MACHINE CONSTANTS. */
/*     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18. */
/*     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. */
/*     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND */
/*     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR */
/*     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. */
/*     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. */
/*     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). */
/*     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU. */
/* ----------------------------------------------------------------------- */
/* Computing MAX */
    r__1 = r1mach_(&c__4);
    tol = dmax(r__1,(float)1e-18);
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
    fnul = (dig - (float)3.) * (float)6. + (float)10.;
    az = c_abs(z__);
/* ----------------------------------------------------------------------- */
/*     TEST FOR RANGE */
/* ----------------------------------------------------------------------- */
    aa = (float).5 / tol;
    bb = i1mach_(&c__9) * (float).5;
    aa = dmin(aa,bb);
    if (az > aa) {
	goto L140;
    }
    fn = *fnu + (*n - 1);
    if (fn > aa) {
	goto L140;
    }
    aa = sqrt(aa);
    if (az > aa) {
	*ierr = 3;
    }
    if (fn > aa) {
	*ierr = 3;
    }
    zn.r = z__->r, zn.i = z__->i;
    csgn.r = cone.r, csgn.i = cone.i;
    if (xx >= (float)0.) {
	goto L40;
    }
    q__1.r = -z__->r, q__1.i = -z__->i;
    zn.r = q__1.r, zn.i = q__1.i;
/* ----------------------------------------------------------------------- */
/*     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE */
/*     WHEN FNU IS LARGE */
/* ----------------------------------------------------------------------- */
    inu = *fnu;
    arg = (*fnu - inu) * pi;
    if (yy < (float)0.) {
	arg = -arg;
    }
    s1 = cos(arg);
    s2 = sin(arg);
    q__1.r = s1, q__1.i = s2;
    csgn.r = q__1.r, csgn.i = q__1.i;
    if (inu % 2 == 1) {
	q__1.r = -csgn.r, q__1.i = -csgn.i;
	csgn.r = q__1.r, csgn.i = q__1.i;
    }
L40:
    cbinu_(&zn, fnu, kode, n, &cy[1], nz, &rl, &fnul, &tol, &elim, &alim);
    if (*nz < 0) {
	goto L120;
    }
    if (xx >= (float)0.) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/*     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE */
/* ----------------------------------------------------------------------- */
    nn = *n - *nz;
    if (nn == 0) {
	return 0;
    }
    rtol = (float)1. / tol;
    ascle = r1mach_(&c__1) * rtol * (float)1e3;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       CY(I) = CY(I)*CSGN */
	i__2 = i__;
	zn.r = cy[i__2].r, zn.i = cy[i__2].i;
	aa = zn.r;
	bb = r_imag(&zn);
	atol = (float)1.;
/* Computing MAX */
	r__1 = dabs(aa), r__2 = dabs(bb);
	if (dmax(r__1,r__2) > ascle) {
	    goto L55;
	}
	q__2.r = rtol, q__2.i = (float)0.;
	q__1.r = zn.r * q__2.r - zn.i * q__2.i, q__1.i = zn.r * q__2.i + zn.i 
		* q__2.r;
	zn.r = q__1.r, zn.i = q__1.i;
	atol = tol;
L55:
	q__1.r = zn.r * csgn.r - zn.i * csgn.i, q__1.i = zn.r * csgn.i + zn.i 
		* csgn.r;
	zn.r = q__1.r, zn.i = q__1.i;
	i__2 = i__;
	q__2.r = atol, q__2.i = (float)0.;
	q__1.r = zn.r * q__2.r - zn.i * q__2.i, q__1.i = zn.r * q__2.i + zn.i 
		* q__2.r;
	cy[i__2].r = q__1.r, cy[i__2].i = q__1.i;
	q__1.r = -csgn.r, q__1.i = -csgn.i;
	csgn.r = q__1.r, csgn.i = q__1.i;
/* L50: */
    }
    return 0;
L120:
    if (*nz == -2) {
	goto L130;
    }
    *nz = 0;
    *ierr = 2;
    return 0;
L130:
    *nz = 0;
    *ierr = 5;
    return 0;
L140:
    *nz = 0;
    *ierr = 4;
    return 0;
} /* cbesi_ */

