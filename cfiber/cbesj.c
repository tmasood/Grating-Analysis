/* cbesj.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/
#include <stdio.h>
#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__5 = 5;
static integer c__11 = 11;
static integer c__9 = 9;
static integer c__1 = 1;

/* DECK CBESJ */
/* Subroutine */ int cbesj(complex *z__, real *fnu, integer *kode,
			    integer *n, complex *cy, integer *nz,
			    integer *ierr)
{
    /* Initialized data */

    static real hpi = (float)1.57079632679489662;

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    complex q__1, q__2;

    /* Builtin functions */
    double r_imag(), c_abs(), sqrt(), cos(), sin();
    void r_cnjg();

    /* Local variables */
    static real alim, elim;
    static complex csgn;
    static real atol;
    static integer inuh;
    static real fnul, rtol;
    static integer i__, k;
    static real ascle;
    extern /* Subroutine */ int cbinu_(complex *, real *, integer *,
				       integer *, complex *, integer *,
				       real *, real *, real *, real *,
				       real *);
    static integer k1, k2;
    extern integer i1mach_();
    static real r1, r2;
    extern doublereal r1mach_();
    static real aa, bb;
    static complex ci;
    static real fn;
    static integer nl;
    static real az;
    static integer ir;
    static real rl;
    static complex zn;
    static real yy, dig, arg, r1m5;
    static integer inu;
    static real tol;

/* ***BEGIN PROLOGUE  CBESJ */
/* ***PURPOSE  Compute a sequence of the Bessel functions J(a,z) for */
/*            complex argument z and real nonnegative orders a=b,b+1, */
/*            b+2,... where b>0.  A scaling option is available to */
/*            help avoid overflow. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10A4 */
/* ***TYPE      COMPLEX (CBESJ-C, ZBESJ-C) */
/* ***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT, */
/*             BESSEL FUNCTIONS OF THE FIRST KIND, J BESSEL FUNCTIONS */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*         On KODE=1, CBESJ computes an N member sequence of complex */
/*         Bessel functions CY(L)=J(FNU+L-1,Z) for real nonnegative */
/*         orders FNU+L-1, L=1,...,N and complex Z in the cut plane */
/*         -pi<arg(Z)<=pi.  On KODE=2, CBESJ returns the scaled functions */

/*            CY(L) = exp(-abs(Y))*J(FNU+L-1,Z),  L=1,...,N and Y=Im(Z) */

/*         which remove the exponential growth in both the upper and */
/*         lower half planes as Z goes to infinity.  Definitions and */
/*         notation are found in the NBS Handbook of Mathematical */
/*         Functions (Ref. 1). */

/*         Input */
/*           Z      - Argument of type COMPLEX */
/*           FNU    - Initial order of type REAL, FNU>=0 */
/*           KODE   - A parameter to indicate the scaling option */
/*                    KODE=1  returns */
/*                            CY(L)=J(FNU+L-1,Z), L=1,...,N */
/*                        =2  returns */
/*                            CY(L)=J(FNU+L-1,Z)*exp(-abs(Y)), L=1,...,N */
/*                            where Y=Im(Z) */
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
/*                            (Im(Z) too large on KODE=1) */
/*                    IERR=3  Precision warning - COMPUTATION COMPLETED */
/*                            (Result has half precision or less */
/*                            because abs(Z) or FNU+N-1 is large) */
/*                    IERR=4  Precision error   - NO COMPUTATION */
/*                            (Result has no precision because */
/*                            abs(Z) or FNU+N-1 is too large) */
/*                    IERR=5  Algorithmic error - NO COMPUTATION */
/*                            (Termination condition not met) */

/* *Long Description: */

/*         The computation is carried out by the formulae */

/*            J(a,z) = exp( a*pi*i/2)*I(a,-i*z),  Im(z)>=0 */

/*            J(a,z) = exp(-a*pi*i/2)*I(a, i*z),  Im(z)<0 */

/*         where the I Bessel function is computed as described in the */
/*         prologue to CBESI. */

/*         For negative orders, the formula */

/*            J(-a,z) = J(a,z)*cos(a*pi) - Y(a,z)*sin(a*pi) */

/*         can be used.  However, for large orders close to integers, the */
/*         the function changes radically.  When a is a large positive */
/*         integer, the magnitude of J(-a,z)=J(a,z)*cos(a*pi) is a */
/*         large negative power of ten.  But when a is not an integer, */
/*         Y(a,z) dominates in magnitude with a large positive power of */
/*         ten and the most that the second term can be reduced is by */
/*         unit roundoff from the coefficient.  Thus, wide changes can */
/*         occur within unit roundoff of a large integer for a.  Here, */
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
/* ***END PROLOGUE  CBESJ */

    /* Parameter adjustments */
    --cy;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  CBESJ */
    *ierr = 0;
    *nz = 0;
    if (*fnu < (float)0.)
      {
	*ierr = 1;
      }

    if (*kode < 1 || *kode > 2)
      {
	*ierr = 1;
      }

    if (*n < 1)
      {
	*ierr = 1;
      }

    if (*ierr != 0)
      {
	return 0;
      }
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
    ci.r = (float)0., ci.i = (float)1.;
    yy = r_imag(z__);
    az = c_abs(z__);
/* ----------------------------------------------------------------------- */
/*     TEST FOR RANGE */
/* ----------------------------------------------------------------------- */
    aa = (float).5 / tol;
    bb = i1mach_(&c__9) * (float).5;
    aa = dmin(aa,bb);
    fn = *fnu + (*n - 1);
    if (az > aa) {
	goto L140;
    }
    if (fn > aa) 
      {
	goto L140;
      }
    aa = sqrt(aa);
    if (az > aa)
      {
	*ierr = 3;
      }
    if (fn > aa)
      {
	*ierr = 3;
      }
/* ----------------------------------------------------------------------- */
/*     CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE */
/*     WHEN FNU IS LARGE */
/* ----------------------------------------------------------------------- */
    inu = *fnu;
    inuh = inu / 2;
    ir = inu - (inuh << 1);
    arg = (*fnu - (inu - ir)) * hpi;
    r1 = cos(arg);
    r2 = sin(arg);
    q__1.r = r1, q__1.i = r2;
    csgn.r = q__1.r, csgn.i = q__1.i;
    if (inuh % 2 == 1) {
	q__1.r = -csgn.r, q__1.i = -csgn.i;
	csgn.r = q__1.r, csgn.i = q__1.i;
    }
/* ----------------------------------------------------------------------- */
/*     ZN IS IN THE RIGHT HALF PLANE */
/* ----------------------------------------------------------------------- */
    q__2.r = -z__->r, q__2.i = -z__->i;

    zn.r = (q__2.r * ci.r) - (q__2.i * ci.i);
    zn.i = (q__2.r * ci.i) + (q__2.i * ci.r);

    if (yy >= (float)0.)
      {
	goto L40;
      }

    q__1.r = -zn.r, q__1.i = -zn.i;
    zn.r = q__1.r, zn.i = q__1.i;
    r_cnjg(&q__1, &csgn);
    csgn.r = q__1.r, csgn.i = q__1.i;
    r_cnjg(&q__1, &ci);
    ci.r = q__1.r, ci.i = q__1.i;
L40:
    cbinu_(&zn, fnu, kode, n, &cy[1], nz, &rl, &fnul, &tol, &elim, &alim);
    if (*nz < 0) {
	goto L120;
    }
    nl = *n - *nz;
    if (nl == 0) {
	return 0;
    }
    rtol = (float)1. / tol;
    ascle = r1mach_(&c__1) * rtol * (float)1e3;
    i__1 = nl;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       CY(I)=CY(I)*CSGN */
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
	q__1.r = csgn.r * ci.r - csgn.i * ci.i, q__1.i = csgn.r * ci.i + 
		csgn.i * ci.r;
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
} /* cbesj_ */

