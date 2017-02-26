/* cbesk.f -- translated by f2c (version 20000121).
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
static integer c__2 = 2;

/* DECK CBESK */
/* Subroutine */ int cbesk(complex *z__, real *fnu, integer *kode,
			    integer *n, complex *cy, integer *nz,
			    integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Builtin functions */
    double r_imag(), c_abs(), sqrt(), log();

    /* Local variables */
    static real alim, elim, fnul;
    static integer k;
    extern /* Subroutine */ int cacon_(complex *, real *, integer *,
			    integer *, integer *, complex *,
			    integer *, real *, real *, real *, real *,
				       real *);
    extern int cbknu_(complex *, real *, integer *,
			    integer *, complex *, integer *,
			    real *, real *, real *);
    extern int cbunk_(complex *, real *, integer *,
			    integer *, integer *, complex *,
			    integer *, real *, real *,
			    real *);
    extern int cuoik_(complex *, real *, integer *, integer *,
		      integer *, complex *, integer *, real *,
		      real *, real *);
    static integer k1, k2;
    extern integer i1mach_();
    extern doublereal r1mach_();
    static real aa, bb, fn, az;
    static integer nn;
    static real rl;
    static integer mr, nw;
    static real xx, yy, dig, arg, aln, r1m5, ufl;
    static integer nuf;
    static real tol;

/* ***BEGIN PROLOGUE  CBESK */
/* ***PURPOSE  Compute a sequence of the Bessel functions K(a,z) for */
/*            complex argument z and real nonnegative orders a=b,b+1, */
/*            b+2,... where b>0.  A scaling option is available to */
/*            help avoid overflow. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10B4 */
/* ***TYPE      COMPLEX (CBESK-C, ZBESK-C) */
/* ***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT, K BESSEL FUNCTIONS, */
/*             MODIFIED BESSEL FUNCTIONS */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*         On KODE=1, CBESK computes an N member sequence of complex */
/*         Bessel functions CY(L)=K(FNU+L-1,Z) for real nonnegative */
/*         orders FNU+L-1, L=1,...,N and complex Z.NE.0 in the cut */
/*         plane -pi<arg(Z)<=pi.  On KODE=2, CBESJ returns the scaled */
/*         functions */

/*            CY(L) = exp(Z)*K(FNU+L-1,Z),  L=1,...,N */

/*         which remove the exponential growth in both the left and */
/*         right half planes as Z goes to infinity.  Definitions and */
/*         notation are found in the NBS Handbook of Mathematical */
/*         Functions (Ref. 1). */

/*         Input */
/*           Z      - Nonzero argument of type COMPLEX */
/*           FNU    - Initial order of type REAL, FNU>=0 */
/*           KODE   - A parameter to indicate the scaling option */
/*                    KODE=1  returns */
/*                            CY(L)=K(FNU+L-1,Z), L=1,...,N */
/*                        =2  returns */
/*                            CY(L)=K(FNU+L-1,Z)*EXP(Z), L=1,...,N */
/*           N      - Number of terms in the sequence, N>=1 */

/*         Output */
/*           CY     - Result vector of type COMPLEX */
/*           NZ     - Number of underflows set to zero */
/*                    NZ=0    Normal return */
/*                    NZ>0    CY(L)=0 for NZ values of L (if Re(Z)>0 */
/*                            then CY(L)=0 for L=1,...,NZ; in the */
/*                            complementary half plane the underflows */
/*                            may not be in an uninterrupted sequence) */
/*           IERR   - Error flag */
/*                    IERR=0  Normal return     - COMPUTATION COMPLETED */
/*                    IERR=1  Input error       - NO COMPUTATION */
/*                    IERR=2  Overflow          - NO COMPUTATION */
/*                            (abs(Z) too small and/or FNU+N-1 */
/*                            too large) */
/*                    IERR=3  Precision warning - COMPUTATION COMPLETED */
/*                            (Result has half precision or less */
/*                            because abs(Z) or FNU+N-1 is large) */
/*                    IERR=4  Precision error   - NO COMPUTATION */
/*                            (Result has no precision because */
/*                            abs(Z) or FNU+N-1 is too large) */
/*                    IERR=5  Algorithmic error - NO COMPUTATION */
/*                            (Termination condition not met) */

/* *Long Description: */

/*         Equations of the reference are implemented to compute K(a,z) */
/*         for small orders a and a+1 in the right half plane Re(z)>=0. */
/*         Forward recurrence generates higher orders.  The formula */

/*            K(a,z*exp((t)) = exp(-t)*K(a,z) - t*I(a,z),  Re(z)>0 */
/*                         t = i*pi or -i*pi */

/*         continues K to the left half plane. */

/*         For large orders, K(a,z) is computed by means of its uniform */
/*         asymptotic expansion. */

/*         For negative orders, the formula */

/*            K(-a,z) = K(a,z) */

/*         can be used. */

/*         CBESK assumes that a significant digit sinh function is */
/*         available. */

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

/* ***ROUTINES CALLED  CACON, CBKNU, CBUNK, CUOIK, I1MACH, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   890801  REVISION DATE from Version 3.2 */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   920128  Category corrected.  (WRB) */
/*   920811  Prologue revised.  (DWL) */
/* ***END PROLOGUE  CBESK */

/* ***FIRST EXECUTABLE STATEMENT  CBESK */
    /* Parameter adjustments */
    --cy;

    /* Function Body */
    *ierr = 0;
    *nz = 0;
    xx = z__->r;
    yy = r_imag(z__);
    if (yy == (float)0. && xx == (float)0.) {
	*ierr = 1;
    }
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
    nn = *n;
/* ----------------------------------------------------------------------- */
/*     SET PARAMETERS RELATED TO MACHINE CONSTANTS. */
/*     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18. */
/*     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT. */
/*     EXP(-ELIM).LT.EXP(-ALIM)=EXP(-ELIM)/TOL    AND */
/*     EXP(ELIM).GT.EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR */
/*     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE. */
/*     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z. */
/*     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG). */
/*     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU */
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
    fnul = (dig - (float)3.) * (float)6. + (float)10.;
    rl = dig * (float)1.2 + (float)3.;
    az = c_abs(z__);
    fn = *fnu + (nn - 1);
/* ----------------------------------------------------------------------- */
/*     TEST FOR RANGE */
/* ----------------------------------------------------------------------- */
    aa = (float).5 / tol;
    bb = i1mach_(&c__9) * (float).5;
    aa = dmin(aa,bb);
    if (az > aa) {
	goto L210;
    }
    if (fn > aa) {
	goto L210;
    }
    aa = sqrt(aa);
    if (az > aa) {
	*ierr = 3;
    }
    if (fn > aa) {
	*ierr = 3;
    }
/* ----------------------------------------------------------------------- */
/*     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE */
/* ----------------------------------------------------------------------- */
/*     UFL = EXP(-ELIM) */
    ufl = r1mach_(&c__1) * (float)1e3;
    if (az < ufl) {
	goto L180;
    }
    if (*fnu > fnul) {
	goto L80;
    }
    if (fn <= (float)1.) {
	goto L60;
    }
    if (fn > (float)2.) {
	goto L50;
    }
    if (az > tol) {
	goto L60;
    }
    arg = az * (float).5;
    aln = -fn * log(arg);
    if (aln > elim) {
	goto L180;
    }
    goto L60;
L50:
    cuoik_(z__, fnu, kode, &c__2, &nn, &cy[1], &nuf, &tol, &elim, &alim);
    if (nuf < 0) {
	goto L180;
    }
    *nz += nuf;
    nn -= nuf;
/* ----------------------------------------------------------------------- */
/*     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK */
/*     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I */
/* ----------------------------------------------------------------------- */
    if (nn == 0) {
	goto L100;
    }
L60:
    if (xx < (float)0.) {
	goto L70;
    }
/* ----------------------------------------------------------------------- */
/*     RIGHT HALF PLANE COMPUTATION, REAL(Z).GE.0. */
/* ----------------------------------------------------------------------- */
    cbknu_(z__, fnu, kode, &nn, &cy[1], &nw, &tol, &elim, &alim);
    if (nw < 0) {
	goto L200;
    }
    *nz = nw;
    return 0;
/* ----------------------------------------------------------------------- */
/*     LEFT HALF PLANE COMPUTATION */
/*     PI/2.LT.ARG(Z).LE.PI AND -PI.LT.ARG(Z).LT.-PI/2. */
/* ----------------------------------------------------------------------- */
L70:
    if (*nz != 0) {
	goto L180;
    }
    mr = 1;
    if (yy < (float)0.) {
	mr = -1;
    }
    cacon_(z__, fnu, kode, &mr, &nn, &cy[1], &nw, &rl, &fnul, &tol, &elim, &
	    alim);
    if (nw < 0) {
	goto L200;
    }
    *nz = nw;
    return 0;
/* ----------------------------------------------------------------------- */
/*     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU.GT.FNUL */
/* ----------------------------------------------------------------------- */
L80:
    mr = 0;
    if (xx >= (float)0.) {
	goto L90;
    }
    mr = 1;
    if (yy < (float)0.) {
	mr = -1;
    }
L90:
    cbunk_(z__, fnu, kode, &mr, &nn, &cy[1], &nw, &tol, &elim, &alim);
    if (nw < 0) {
	goto L200;
    }
    *nz += nw;
    return 0;
L100:
    if (xx < (float)0.) {
	goto L180;
    }
    return 0;
L180:
    *nz = 0;
    *ierr = 2;
    return 0;
L200:
    if (nw == -1) {
	goto L180;
    }
    *nz = 0;
    *ierr = 5;
    return 0;
L210:
    *nz = 0;
    *ierr = 4;
    return 0;
} /* cbesk_ */

