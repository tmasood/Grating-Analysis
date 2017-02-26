/* cbesy.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__4 = 4;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__5 = 5;

/* DECK CBESY */
/* Subroutine */ int cbesy(complex *z__, real *fnu, integer *kode,
			   integer *n, complex *cy, integer *nz,
			   complex *cwrk,  integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    real r__1, r__2;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    double r_imag(), cos(), sin(), exp();
    void r_cnjg();

    /* Local variables */
    static real elim, atol, rtol;
    static integer i__, k;
    extern /* Subroutine */ int cbesh_(complex *, real *, integer *,
			    integer *, integer *, complex *, 
			    integer *, integer *);
    static real ascle;
    static complex c1, c2;
    static integer k1, k2;
    extern integer i1mach_();
    static real r1, r2;
    extern doublereal r1mach_();
    static real aa, bb;
    static complex ex;
    static real ey;
    static complex zu, zv;
    static real xx, yy;
    static integer nz1, nz2;
    static complex hci;
    static real r1m5, tay, tol;

/* ***BEGIN PROLOGUE  CBESY */
/* ***PURPOSE  Compute a sequence of the Bessel functions Y(a,z) for */
/*            complex argument z and real nonnegative orders a=b,b+1, */
/*            b+2,... where b>0.  A scaling option is available to */
/*            help avoid overflow. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C10A4 */
/* ***TYPE      COMPLEX (CBESY-C, ZBESY-C) */
/* ***KEYWORDS  BESSEL FUNCTIONS OF COMPLEX ARGUMENT, */
/*             BESSEL FUNCTIONS OF SECOND KIND, WEBER'S FUNCTION, */
/*             Y BESSEL FUNCTIONS */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*         On KODE=1, CBESY computes an N member sequence of complex */
/*         Bessel functions CY(L)=Y(FNU+L-1,Z) for real nonnegative */
/*         orders FNU+L-1, L=1,...,N and complex Z in the cut plane */
/*         -pi<arg(Z)<=pi.  On KODE=2, CBESY returns the scaled */
/*         functions */

/*            CY(L) = exp(-abs(Y))*Y(FNU+L-1,Z),  L=1,...,N, Y=Im(Z) */

/*         which remove the exponential growth in both the upper and */
/*         lower half planes as Z goes to infinity.  Definitions and */
/*         notation are found in the NBS Handbook of Mathematical */
/*         Functions (Ref. 1). */

/*         Input */
/*           Z      - Nonzero argument of type COMPLEX */
/*           FNU    - Initial order of type REAL, FNU>=0 */
/*           KODE   - A parameter to indicate the scaling option */
/*                    KODE=1  returns */
/*                            CY(L)=Y(FNU+L-1,Z), L=1,...,N */
/*                        =2  returns */
/*                            CY(L)=Y(FNU+L-1,Z)*exp(-abs(Y)), L=1,...,N */
/*                            where Y=Im(Z) */
/*           N      - Number of terms in the sequence, N>=1 */
/*           CWRK   - A work vector of type COMPLEX and dimension N */

/*         Output */
/*           CY     - Result vector of type COMPLEX */
/*           NZ     - Number of underflows set to zero */
/*                    NZ=0    Normal return */
/*                    NZ>0    CY(L)=0 for NZ values of L, usually on */
/*                            KODE=2 (the underflows may not be in an */
/*                            uninterrupted sequence) */
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

/*         The computation is carried out by the formula */

/*            Y(a,z) = (H(1,a,z) - H(2,a,z))/(2*i) */

/*         where the Hankel functions are computed as described in CBESH. */

/*         For negative orders, the formula */

/*            Y(-a,z) = Y(a,z)*cos(a*pi) + J(a,z)*sin(a*pi) */

/*         can be used.  However, for large orders close to half odd */
/*         integers the function changes radically.  When a is a large */
/*         positive half odd integer, the magnitude of Y(-a,z)=J(a,z)* */
/*         sin(a*pi) is a large negative power of ten.  But when a is */
/*         not a half odd integer, Y(a,z) dominates in magnitude with a */
/*         large positive power of ten and the most that the second term */
/*         can be reduced is by unit roundoff from the coefficient. */
/*         Thus,  wide changes can occur within unit roundoff of a large */
/*         half odd integer.  Here, large means a>abs(z). */

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

/* ***ROUTINES CALLED  CBESH, I1MACH, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   890801  REVISION DATE from Version 3.2 */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/*   920128  Category corrected.  (WRB) */
/*   920811  Prologue revised.  (DWL) */
/* ***END PROLOGUE  CBESY */

/* ***FIRST EXECUTABLE STATEMENT  CBESY */
    /* Parameter adjustments */
    --cwrk;
    --cy;

    /* Function Body */
    xx = z__->r;
    yy = r_imag(z__);
    *ierr = 0;
    *nz = 0;
    if (xx == (float)0. && yy == (float)0.) {
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
    hci.r = (float)0., hci.i = (float).5;
    cbesh_(z__, fnu, kode, &c__1, n, &cy[1], &nz1, ierr);
    if (*ierr != 0 && *ierr != 3) {
	goto L170;
    }
    cbesh_(z__, fnu, kode, &c__2, n, &cwrk[1], &nz2, ierr);
    if (*ierr != 0 && *ierr != 3) {
	goto L170;
    }
    *nz = min(nz1,nz2);
    if (*kode == 2) {
	goto L60;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	i__4 = i__;
	q__2.r = cwrk[i__3].r - cy[i__4].r, q__2.i = cwrk[i__3].i - cy[i__4]
		.i;
	q__1.r = hci.r * q__2.r - hci.i * q__2.i, q__1.i = hci.r * q__2.i + 
		hci.i * q__2.r;
	cy[i__2].r = q__1.r, cy[i__2].i = q__1.i;
/* L50: */
    }
    return 0;
L60:
/* Computing MAX */
    r__1 = r1mach_(&c__4);
    tol = dmax(r__1,(float)1e-18);
    k1 = i1mach_(&c__12);
    k2 = i1mach_(&c__13);
/* Computing MIN */
    i__1 = abs(k1), i__2 = abs(k2);
    k = min(i__1,i__2);
    r1m5 = r1mach_(&c__5);
/* ----------------------------------------------------------------------- */
/*     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT */
/* ----------------------------------------------------------------------- */
    elim = (k * r1m5 - (float)3.) * (float)2.303;
    r1 = cos(xx);
    r2 = sin(xx);
    q__1.r = r1, q__1.i = r2;
    ex.r = q__1.r, ex.i = q__1.i;
    ey = (float)0.;
    tay = (r__1 = yy + yy, dabs(r__1));
    if (tay < elim) {
	ey = exp(-tay);
    }
    if (yy < (float)0.) {
	goto L90;
    }
    q__2.r = ey, q__2.i = (float)0.;
    q__1.r = ex.r * q__2.r - ex.i * q__2.i, q__1.i = ex.r * q__2.i + ex.i * 
	    q__2.r;
    c1.r = q__1.r, c1.i = q__1.i;
    r_cnjg(&q__1, &ex);
    c2.r = q__1.r, c2.i = q__1.i;
L70:
    *nz = 0;
    rtol = (float)1. / tol;
    ascle = r1mach_(&c__1) * rtol * (float)1e3;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*       CY(I) = HCI*(C2*CWRK(I)-C1*CY(I)) */
	i__2 = i__;
	zv.r = cwrk[i__2].r, zv.i = cwrk[i__2].i;
	aa = zv.r;
	bb = r_imag(&zv);
	atol = (float)1.;
/* Computing MAX */
	r__1 = dabs(aa), r__2 = dabs(bb);
	if (dmax(r__1,r__2) > ascle) {
	    goto L75;
	}
	q__2.r = rtol, q__2.i = (float)0.;
	q__1.r = zv.r * q__2.r - zv.i * q__2.i, q__1.i = zv.r * q__2.i + zv.i 
		* q__2.r;
	zv.r = q__1.r, zv.i = q__1.i;
	atol = tol;
L75:
	q__2.r = zv.r * c2.r - zv.i * c2.i, q__2.i = zv.r * c2.i + zv.i * 
		c2.r;
	q__1.r = q__2.r * hci.r - q__2.i * hci.i, q__1.i = q__2.r * hci.i + 
		q__2.i * hci.r;
	zv.r = q__1.r, zv.i = q__1.i;
	q__2.r = atol, q__2.i = (float)0.;
	q__1.r = zv.r * q__2.r - zv.i * q__2.i, q__1.i = zv.r * q__2.i + zv.i 
		* q__2.r;
	zv.r = q__1.r, zv.i = q__1.i;
	i__2 = i__;
	zu.r = cy[i__2].r, zu.i = cy[i__2].i;
	aa = zu.r;
	bb = r_imag(&zu);
	atol = (float)1.;
/* Computing MAX */
	r__1 = dabs(aa), r__2 = dabs(bb);
	if (dmax(r__1,r__2) > ascle) {
	    goto L85;
	}
	q__2.r = rtol, q__2.i = (float)0.;
	q__1.r = zu.r * q__2.r - zu.i * q__2.i, q__1.i = zu.r * q__2.i + zu.i 
		* q__2.r;
	zu.r = q__1.r, zu.i = q__1.i;
	atol = tol;
L85:
	q__2.r = zu.r * c1.r - zu.i * c1.i, q__2.i = zu.r * c1.i + zu.i * 
		c1.r;
	q__1.r = q__2.r * hci.r - q__2.i * hci.i, q__1.i = q__2.r * hci.i + 
		q__2.i * hci.r;
	zu.r = q__1.r, zu.i = q__1.i;
	q__2.r = atol, q__2.i = (float)0.;
	q__1.r = zu.r * q__2.r - zu.i * q__2.i, q__1.i = zu.r * q__2.i + zu.i 
		* q__2.r;
	zu.r = q__1.r, zu.i = q__1.i;
	i__2 = i__;
	q__1.r = zv.r - zu.r, q__1.i = zv.i - zu.i;
	cy[i__2].r = q__1.r, cy[i__2].i = q__1.i;
	i__2 = i__;
	if (cy[i__2].r == (float)0. && cy[i__2].i == (float)0. && ey == (
		float)0.) {
	    ++(*nz);
	}
/* L80: */
    }
    return 0;
L90:
    c1.r = ex.r, c1.i = ex.i;
    r_cnjg(&q__2, &ex);
    q__3.r = ey, q__3.i = (float)0.;
    q__1.r = q__2.r * q__3.r - q__2.i * q__3.i, q__1.i = q__2.r * q__3.i + 
	    q__2.i * q__3.r;
    c2.r = q__1.r, c2.i = q__1.i;
    goto L70;
L170:
    *nz = 0;
    return 0;
} /* cbesy_ */

