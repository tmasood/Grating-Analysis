#include "f2c.h"

/* Subroutine */ int zlahqr_(logical *wantt, logical *wantz, integer *n, 
	integer *ilo, integer *ihi, doublecomplex *h, integer *ldh, 
	doublecomplex *w, integer *iloz, integer *ihiz, doublecomplex *z, 
	integer *ldz, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    ZLAHQR is an auxiliary routine called by ZHSEQR to update the   
    eigenvalues and Schur decomposition already computed by ZHSEQR, by   
    dealing with the Hessenberg submatrix in rows and columns ILO to IHI. 
  

    Arguments   
    =========   

    WANTT   (input) LOGICAL   
            = .TRUE. : the full Schur form T is required;   
            = .FALSE.: only eigenvalues are required.   

    WANTZ   (input) LOGICAL   
            = .TRUE. : the matrix of Schur vectors Z is required;   
            = .FALSE.: Schur vectors are not required.   

    N       (input) INTEGER   
            The order of the matrix H.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            It is assumed that H is already upper triangular in rows and 
  
            columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1). 
  
            ZLAHQR works primarily with the Hessenberg submatrix in rows 
  
            and columns ILO to IHI, but applies transformations to all of 
  
            H if WANTT is .TRUE..   
            1 <= ILO <= max(1,IHI); IHI <= N.   

    H       (input/output) COMPLEX*16 array, dimension (LDH,N)   
            On entry, the upper Hessenberg matrix H.   
            On exit, if WANTT is .TRUE., H is upper triangular in rows   
            and columns ILO:IHI, with any 2-by-2 diagonal blocks in   
            standard form. If WANTT is .FALSE., the contents of H are   
            unspecified on exit.   

    LDH     (input) INTEGER   
            The leading dimension of the array H. LDH >= max(1,N).   

    W       (output) COMPLEX*16 array, dimension (N)   
            The computed eigenvalues ILO to IHI are stored in the   
            corresponding elements of W. If WANTT is .TRUE., the   
            eigenvalues are stored in the same order as on the diagonal   
            of the Schur form returned in H, with W(i) = H(i,i).   

    ILOZ    (input) INTEGER   
    IHIZ    (input) INTEGER   
            Specify the rows of Z to which transformations must be   
            applied if WANTZ is .TRUE..   
            1 <= ILOZ <= ILO; IHI <= IHIZ <= N.   

    Z       (input/output) COMPLEX*16 array, dimension (LDZ,N)   
            If WANTZ is .TRUE., on entry Z must contain the current   
            matrix Z of transformations accumulated by ZHSEQR, and on   
            exit Z has been updated; transformations are applied only to 
  
            the submatrix Z(ILOZ:IHIZ,ILO:IHI).   
            If WANTZ is .FALSE., Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z. LDZ >= max(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            > 0: if INFO = i, ZLAHQR failed to compute all the   
                 eigenvalues ILO to IHI in a total of 30*(IHI-ILO+1)   
                 iterations; elements i+1:ihi of W contain those   
                 eigenvalues which have been successfully computed.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__2 = 2;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer h_dim1, h_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    doublecomplex z__1, z__2, z__3, z__4;
    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_sqrt(doublecomplex *, doublecomplex *), d_cnjg(doublecomplex *, 
	    doublecomplex *);
    /* Local variables */
    static doublereal unfl, ovfl;
    static doublecomplex temp;
    static integer i, j, k, l, m;
    static doublereal s;
    static doublecomplex t, u, v[2], x, y;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static doublereal rtemp;
    static integer i1, i2;
    static doublereal rwork[1];
    static doublecomplex t1;
    static doublereal t2;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static doublecomplex v2;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    static doublereal h10;
    static doublecomplex h11;
    static doublereal h21;
    static doublecomplex h22;
    static integer nh;
    extern doublereal dlamch_(char *);
    static integer nz;
    extern /* Subroutine */ int zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *);
    extern /* Double Complex */ VOID zladiv_(doublecomplex *, doublecomplex *,
	     doublecomplex *);
    extern doublereal zlanhs_(char *, integer *, doublecomplex *, integer *, 
	    doublereal *);
    static doublereal smlnum;
    static doublecomplex h11s;
    static integer itn, its;
    static doublereal ulp;
    static doublecomplex sum;
    static doublereal tst1;



#define V(I) v[(I)]
#define RWORK(I) rwork[(I)]
#define W(I) w[(I)-1]

#define H(I,J) h[(I)-1 + ((J)-1)* ( *ldh)]
#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    *info = 0;

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }
    if (*ilo == *ihi) {
	i__1 = *ilo;
	i__2 = *ilo + *ilo * h_dim1;
	W(*ilo).r = H(*ilo,*ilo).r, W(*ilo).i = H(*ilo,*ilo).i;
	return 0;
    }

    nh = *ihi - *ilo + 1;
    nz = *ihiz - *iloz + 1;

/*     Set machine-dependent constants for the stopping criterion.   
       If norm(H) <= sqrt(OVFL), overflow should not occur. */

    unfl = dlamch_("Safe minimum");
    ovfl = 1. / unfl;
    dlabad_(&unfl, &ovfl);
    ulp = dlamch_("Precision");
    smlnum = unfl * (nh / ulp);

/*     I1 and I2 are the indices of the first row and last column of H   
       to which transformations must be applied. If eigenvalues only are 
  
       being computed, I1 and I2 are set inside the main loop. */

    if (*wantt) {
	i1 = 1;
	i2 = *n;
    }

/*     ITN is the total number of QR iterations allowed. */

    itn = nh * 30;

/*     The main loop begins here. I is the loop index and decreases from 
  
       IHI to ILO in steps of 1. Each iteration of the loop works   
       with the active submatrix in rows and columns L to I.   
       Eigenvalues I+1 to IHI have already converged. Either L = ILO, or 
  
       H(L,L-1) is negligible so that the matrix splits. */

    i = *ihi;
L10:
    if (i < *ilo) {
	goto L130;
    }

/*     Perform QR iterations on rows and columns ILO to I until a   
       submatrix of order 1 splits off at the bottom because a   
       subdiagonal element has become negligible. */

    l = *ilo;
    i__1 = itn;
    for (its = 0; its <= itn; ++its) {

/*        Look for a single small subdiagonal element. */

	i__2 = l + 1;
	for (k = i; k >= l+1; --k) {
	    i__3 = k - 1 + (k - 1) * h_dim1;
	    i__4 = k + k * h_dim1;
	    tst1 = (d__1 = H(k-1,k-1).r, abs(d__1)) + (d__2 = d_imag(&H(k-1,k-1)), abs(d__2)) + ((d__3 = H(k,k).r, abs(
		    d__3)) + (d__4 = d_imag(&H(k,k)), abs(d__4)));
	    if (tst1 == 0.) {
		i__3 = i - l + 1;
		tst1 = zlanhs_("1", &i__3, &H(l,l), ldh, rwork)
			;
	    }
	    i__3 = k + (k - 1) * h_dim1;
/* Computing MAX */
	    d__2 = ulp * tst1;
	    if ((d__1 = H(k,k-1).r, abs(d__1)) <= max(d__2,smlnum)) {
		goto L30;
	    }
/* L20: */
	}
L30:
	l = k;
	if (l > *ilo) {

/*           H(L,L-1) is negligible */

	    i__2 = l + (l - 1) * h_dim1;
	    H(l,l-1).r = 0., H(l,l-1).i = 0.;
	}

/*        Exit from loop if a submatrix of order 1 has split off. */

	if (l >= i) {
	    goto L120;
	}

/*        Now the active submatrix is in rows and columns L to I. If 
  
          eigenvalues only are being computed, only the active submatr
ix   
          need be transformed. */

	if (! (*wantt)) {
	    i1 = l;
	    i2 = i;
	}

	if (its == 10 || its == 20) {

/*           Exceptional shift. */

	    i__2 = i + (i - 1) * h_dim1;
	    i__3 = i - 1 + (i - 2) * h_dim1;
	    d__3 = (d__1 = H(i,i-1).r, abs(d__1)) + (d__2 = H(i-1,i-2).r, abs(
		    d__2));
	    t.r = d__3, t.i = 0.;
	} else {

/*           Wilkinson's shift. */

	    i__2 = i + i * h_dim1;
	    t.r = H(i,i).r, t.i = H(i,i).i;
	    i__2 = i - 1 + i * h_dim1;
	    i__3 = i + (i - 1) * h_dim1;
	    d__1 = H(i,i-1).r;
	    z__1.r = d__1 * H(i-1,i).r, z__1.i = d__1 * H(i-1,i).i;
	    u.r = z__1.r, u.i = z__1.i;
	    if (u.r != 0. || u.i != 0.) {
		i__2 = i - 1 + (i - 1) * h_dim1;
		z__2.r = H(i-1,i-1).r - t.r, z__2.i = H(i-1,i-1).i - t.i;
		z__1.r = z__2.r * .5, z__1.i = z__2.i * .5;
		x.r = z__1.r, x.i = z__1.i;
		z__3.r = x.r * x.r - x.i * x.i, z__3.i = x.r * x.i + x.i * 
			x.r;
		z__2.r = z__3.r + u.r, z__2.i = z__3.i + u.i;
		z_sqrt(&z__1, &z__2);
		y.r = z__1.r, y.i = z__1.i;
		if (x.r * y.r + d_imag(&x) * d_imag(&y) < 0.) {
		    z__1.r = -y.r, z__1.i = -y.i;
		    y.r = z__1.r, y.i = z__1.i;
		}
		z__3.r = x.r + y.r, z__3.i = x.i + y.i;
		zladiv_(&z__2, &u, &z__3);
		z__1.r = t.r - z__2.r, z__1.i = t.i - z__2.i;
		t.r = z__1.r, t.i = z__1.i;
	    }
	}

/*        Look for two consecutive small subdiagonal elements. */

	i__2 = l;
	for (m = i - 1; m >= l; --m) {

/*           Determine the effect of starting the single-shift QR 
  
             iteration at row M, and see if this would make H(M,M-
1)   
             negligible. */

	    i__3 = m + m * h_dim1;
	    h11.r = H(m,m).r, h11.i = H(m,m).i;
	    i__3 = m + 1 + (m + 1) * h_dim1;
	    h22.r = H(m+1,m+1).r, h22.i = H(m+1,m+1).i;
	    z__1.r = h11.r - t.r, z__1.i = h11.i - t.i;
	    h11s.r = z__1.r, h11s.i = z__1.i;
	    i__3 = m + 1 + m * h_dim1;
	    h21 = H(m+1,m).r;
	    s = (d__1 = h11s.r, abs(d__1)) + (d__2 = d_imag(&h11s), abs(d__2))
		     + abs(h21);
	    z__1.r = h11s.r / s, z__1.i = h11s.i / s;
	    h11s.r = z__1.r, h11s.i = z__1.i;
	    h21 /= s;
	    V(0).r = h11s.r, V(0).i = h11s.i;
	    V(1).r = h21, V(1).i = 0.;
	    if (m == l) {
		goto L50;
	    }
	    i__3 = m + (m - 1) * h_dim1;
	    h10 = H(m,m-1).r;
	    tst1 = ((d__1 = h11s.r, abs(d__1)) + (d__2 = d_imag(&h11s), abs(
		    d__2))) * ((d__3 = h11.r, abs(d__3)) + (d__4 = d_imag(&
		    h11), abs(d__4)) + ((d__5 = h22.r, abs(d__5)) + (d__6 = 
		    d_imag(&h22), abs(d__6))));
	    if ((d__1 = h10 * h21, abs(d__1)) <= ulp * tst1) {
		goto L50;
	    }
/* L40: */
	}
L50:

/*        Single-shift QR step */

	i__2 = i - 1;
	for (k = m; k <= i-1; ++k) {

/*           The first iteration of this loop determines a reflect
ion G   
             from the vector V and applies it from left and right 
to H,   
             thus creating a nonzero bulge below the subdiagonal. 
  

             Each subsequent iteration determines a reflection G t
o   
             restore the Hessenberg form in the (K-1)th column, an
d thus   
             chases the bulge one step toward the bottom of the ac
tive   
             submatrix.   

             V(2) is always real before the call to ZLARFG, and he
nce   
             after the call T2 ( = T1*V(2) ) is also real. */

	    if (k > m) {
		zcopy_(&c__2, &H(k,k-1), &c__1, v, &c__1);
	    }
	    zlarfg_(&c__2, v, &V(1), &c__1, &t1);
	    if (k > m) {
		i__3 = k + (k - 1) * h_dim1;
		H(k,k-1).r = V(0).r, H(k,k-1).i = V(0).i;
		i__3 = k + 1 + (k - 1) * h_dim1;
		H(k+1,k-1).r = 0., H(k+1,k-1).i = 0.;
	    }
	    v2.r = V(1).r, v2.i = V(1).i;
	    z__1.r = t1.r * v2.r - t1.i * v2.i, z__1.i = t1.r * v2.i + t1.i * 
		    v2.r;
	    t2 = z__1.r;

/*           Apply G from the left to transform the rows of the ma
trix   
             in columns K to I2. */

	    i__3 = i2;
	    for (j = k; j <= i2; ++j) {
		d_cnjg(&z__3, &t1);
		i__4 = k + j * h_dim1;
		z__2.r = z__3.r * H(k,j).r - z__3.i * H(k,j).i, z__2.i = 
			z__3.r * H(k,j).i + z__3.i * H(k,j).r;
		i__5 = k + 1 + j * h_dim1;
		z__4.r = t2 * H(k+1,j).r, z__4.i = t2 * H(k+1,j).i;
		z__1.r = z__2.r + z__4.r, z__1.i = z__2.i + z__4.i;
		sum.r = z__1.r, sum.i = z__1.i;
		i__4 = k + j * h_dim1;
		i__5 = k + j * h_dim1;
		z__1.r = H(k,j).r - sum.r, z__1.i = H(k,j).i - sum.i;
		H(k,j).r = z__1.r, H(k,j).i = z__1.i;
		i__4 = k + 1 + j * h_dim1;
		i__5 = k + 1 + j * h_dim1;
		z__2.r = sum.r * v2.r - sum.i * v2.i, z__2.i = sum.r * v2.i + 
			sum.i * v2.r;
		z__1.r = H(k+1,j).r - z__2.r, z__1.i = H(k+1,j).i - z__2.i;
		H(k+1,j).r = z__1.r, H(k+1,j).i = z__1.i;
/* L60: */
	    }

/*           Apply G from the right to transform the columns of th
e   
             matrix in rows I1 to min(K+2,I).   

   Computing MIN */
	    i__4 = k + 2;
	    i__3 = min(i__4,i);
	    for (j = i1; j <= min(k+2,i); ++j) {
		i__4 = j + k * h_dim1;
		z__2.r = t1.r * H(j,k).r - t1.i * H(j,k).i, z__2.i = t1.r * 
			H(j,k).i + t1.i * H(j,k).r;
		i__5 = j + (k + 1) * h_dim1;
		z__3.r = t2 * H(j,k+1).r, z__3.i = t2 * H(j,k+1).i;
		z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		sum.r = z__1.r, sum.i = z__1.i;
		i__4 = j + k * h_dim1;
		i__5 = j + k * h_dim1;
		z__1.r = H(j,k).r - sum.r, z__1.i = H(j,k).i - sum.i;
		H(j,k).r = z__1.r, H(j,k).i = z__1.i;
		i__4 = j + (k + 1) * h_dim1;
		i__5 = j + (k + 1) * h_dim1;
		d_cnjg(&z__3, &v2);
		z__2.r = sum.r * z__3.r - sum.i * z__3.i, z__2.i = sum.r * 
			z__3.i + sum.i * z__3.r;
		z__1.r = H(j,k+1).r - z__2.r, z__1.i = H(j,k+1).i - z__2.i;
		H(j,k+1).r = z__1.r, H(j,k+1).i = z__1.i;
/* L70: */
	    }

	    if (*wantz) {

/*              Accumulate transformations in the matrix Z */

		i__3 = *ihiz;
		for (j = *iloz; j <= *ihiz; ++j) {
		    i__4 = j + k * z_dim1;
		    z__2.r = t1.r * Z(j,k).r - t1.i * Z(j,k).i, z__2.i = 
			    t1.r * Z(j,k).i + t1.i * Z(j,k).r;
		    i__5 = j + (k + 1) * z_dim1;
		    z__3.r = t2 * Z(j,k+1).r, z__3.i = t2 * Z(j,k+1).i;
		    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
		    sum.r = z__1.r, sum.i = z__1.i;
		    i__4 = j + k * z_dim1;
		    i__5 = j + k * z_dim1;
		    z__1.r = Z(j,k).r - sum.r, z__1.i = Z(j,k).i - sum.i;
		    Z(j,k).r = z__1.r, Z(j,k).i = z__1.i;
		    i__4 = j + (k + 1) * z_dim1;
		    i__5 = j + (k + 1) * z_dim1;
		    d_cnjg(&z__3, &v2);
		    z__2.r = sum.r * z__3.r - sum.i * z__3.i, z__2.i = sum.r *
			     z__3.i + sum.i * z__3.r;
		    z__1.r = Z(j,k+1).r - z__2.r, z__1.i = Z(j,k+1).i - z__2.i;
		    Z(j,k+1).r = z__1.r, Z(j,k+1).i = z__1.i;
/* L80: */
		}
	    }

	    if (k == m && m > l) {

/*              If the QR step was started at row M > L becaus
e two   
                consecutive small subdiagonals were found, the
n extra   
                scaling must be performed to ensure that H(M,M
-1) remains   
                real. */

		z__1.r = 1. - t1.r, z__1.i = 0. - t1.i;
		temp.r = z__1.r, temp.i = z__1.i;
		d__2 = temp.r;
		d__3 = d_imag(&temp);
		d__1 = dlapy2_(&d__2, &d__3);
		z__1.r = temp.r / d__1, z__1.i = temp.i / d__1;
		temp.r = z__1.r, temp.i = z__1.i;
		i__3 = m + 1 + m * h_dim1;
		i__4 = m + 1 + m * h_dim1;
		d_cnjg(&z__2, &temp);
		z__1.r = H(m+1,m).r * z__2.r - H(m+1,m).i * z__2.i, z__1.i = H(m+1,m).r * z__2.i + H(m+1,m).i * z__2.r;
		H(m+1,m).r = z__1.r, H(m+1,m).i = z__1.i;
		if (m + 2 <= i) {
		    i__3 = m + 2 + (m + 1) * h_dim1;
		    i__4 = m + 2 + (m + 1) * h_dim1;
		    z__1.r = H(m+2,m+1).r * temp.r - H(m+2,m+1).i * temp.i, z__1.i =
			     H(m+2,m+1).r * temp.i + H(m+2,m+1).i * temp.r;
		    H(m+2,m+1).r = z__1.r, H(m+2,m+1).i = z__1.i;
		}
		i__3 = i;
		for (j = m; j <= i; ++j) {
		    if (j != m + 1) {
			if (i2 > j) {
			    i__4 = i2 - j;
			    zscal_(&i__4, &temp, &H(j,j+1), 
				    ldh);
			}
			i__4 = j - i1;
			d_cnjg(&z__1, &temp);
			zscal_(&i__4, &z__1, &H(i1,j), &c__1);
			if (*wantz) {
			    d_cnjg(&z__1, &temp);
			    zscal_(&nz, &z__1, &Z(*iloz,j), &c__1);
			}
		    }
/* L90: */
		}
	    }
/* L100: */
	}

/*        Ensure that H(I,I-1) is real. */

	i__2 = i + (i - 1) * h_dim1;
	temp.r = H(i,i-1).r, temp.i = H(i,i-1).i;
	if (d_imag(&temp) != 0.) {
	    d__1 = temp.r;
	    d__2 = d_imag(&temp);
	    rtemp = dlapy2_(&d__1, &d__2);
	    i__2 = i + (i - 1) * h_dim1;
	    H(i,i-1).r = rtemp, H(i,i-1).i = 0.;
	    z__1.r = temp.r / rtemp, z__1.i = temp.i / rtemp;
	    temp.r = z__1.r, temp.i = z__1.i;
	    if (i2 > i) {
		i__2 = i2 - i;
		d_cnjg(&z__1, &temp);
		zscal_(&i__2, &z__1, &H(i,i+1), ldh);
	    }
	    i__2 = i - i1;
	    zscal_(&i__2, &temp, &H(i1,i), &c__1);
	    if (*wantz) {
		zscal_(&nz, &temp, &Z(*iloz,i), &c__1);
	    }
	}

/* L110: */
    }

/*     Failure to converge in remaining number of iterations */

    *info = i;
    return 0;

L120:

/*     H(I,I-1) is negligible: one eigenvalue has converged. */

    i__1 = i;
    i__2 = i + i * h_dim1;
    W(i).r = H(i,i).r, W(i).i = H(i,i).i;

/*     Decrement number of remaining iterations, and return to start of   
       the main loop with new value of I. */

    itn -= its;
    i = l - 1;
    goto L10;

L130:
    return 0;

/*     End of ZLAHQR */

} /* zlahqr_ */

