#include "f2c.h"

/* Subroutine */ int dsteqr_(char *compz, integer *n, doublereal *d, 
	doublereal *e, doublereal *z, integer *ldz, doublereal *work, integer 
	*info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSTEQR computes all eigenvalues and, optionally, eigenvectors of a   
    symmetric tridiagonal matrix using the implicit QL or QR method.   
    The eigenvectors of a full or band symmetric matrix can also be found 
  
    if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to 
  
    tridiagonal form.   

    Arguments   
    =========   

    COMPZ   (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only.   
            = 'V':  Compute eigenvalues and eigenvectors of the original 
  
                    symmetric matrix.  On entry, Z must contain the   
                    orthogonal matrix used to reduce the original matrix 
  
                    to tridiagonal form.   
            = 'I':  Compute eigenvalues and eigenvectors of the   
                    tridiagonal matrix.  Z is initialized to the identity 
  
                    matrix.   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) DOUBLE PRECISION array, dimension (N)   
            On entry, the diagonal elements of the tridiagonal matrix.   
            On exit, if INFO = 0, the eigenvalues in ascending order.   

    E       (input/output) DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix.   
            On exit, E has been destroyed.   

    Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)   
            On entry, if  COMPZ = 'V', then Z contains the orthogonal   
            matrix used in the reduction to tridiagonal form.   
            On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the   
            orthonormal eigenvectors of the original symmetric matrix,   
            and if COMPZ = 'I', Z contains the orthonormal eigenvectors   
            of the symmetric tridiagonal matrix.   
            If COMPZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            eigenvectors are desired, then  LDZ >= max(1,N).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2)) 
  
            If COMPZ = 'N', then WORK is not referenced.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  the algorithm has failed to find all the eigenvalues in 
  
                  a total of 30*N iterations; if INFO = i, then i   
                  elements of E have not converged to zero; on exit, D   
                  and E contain the elements of a symmetric tridiagonal   
                  matrix which is orthogonally similar to the original   
                  matrix.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static doublereal c_b9 = 0.;
    static doublereal c_b10 = 1.;
    static integer c__0 = 0;
    static integer c__1 = 1;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;
    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Local variables */
    static integer lend, jtot;
    extern /* Subroutine */ int dlae2_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *);
    static doublereal b, c, f, g;
    static integer i, j, k, l, m;
    static doublereal p, r, s;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *);
    static doublereal anorm;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer l1;
    extern /* Subroutine */ int dlaev2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static integer lendm1, lendp1;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static integer ii;
    extern doublereal dlamch_(char *);
    static integer mm, iscale;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *), dlaset_(char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *);
    static doublereal safmin;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal safmax;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    extern doublereal dlanst_(char *, integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int dlasrt_(char *, integer *, doublereal *, 
	    integer *);
    static integer lendsv;
    static doublereal ssfmin;
    static integer nmaxit, icompz;
    static doublereal ssfmax;
    static integer lm1, mm1, nm1;
    static doublereal rt1, rt2, eps;
    static integer lsv;
    static doublereal tst, eps2;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    *info = 0;

    if (lsame_(compz, "N")) {
	icompz = 0;
    } else if (lsame_(compz, "V")) {
	icompz = 1;
    } else if (lsame_(compz, "I")) {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (icompz < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || icompz > 0 && *ldz < max(1,*n)) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSTEQR", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (*n == 1) {
	if (icompz == 2) {
	    Z(1,1) = 1.;
	}
	return 0;
    }

/*     Determine the unit roundoff and over/underflow thresholds. */

    eps = dlamch_("E");
/* Computing 2nd power */
    d__1 = eps;
    eps2 = d__1 * d__1;
    safmin = dlamch_("S");
    safmax = 1. / safmin;
    ssfmax = sqrt(safmax) / 3.;
    ssfmin = sqrt(safmin) / eps2;

/*     Compute the eigenvalues and eigenvectors of the tridiagonal   
       matrix. */

    if (icompz == 2) {
	dlaset_("Full", n, n, &c_b9, &c_b10, &Z(1,1), ldz);
    }

    nmaxit = *n * 30;
    jtot = 0;

/*     Determine where the matrix splits and choose QL or QR iteration   
       for each block, according to whether top or bottom diagonal   
       element is smaller. */

    l1 = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n) {
	goto L160;
    }
    if (l1 > 1) {
	E(l1 - 1) = 0.;
    }
    if (l1 <= nm1) {
	i__1 = nm1;
	for (m = l1; m <= nm1; ++m) {
	    tst = (d__1 = E(m), abs(d__1));
	    if (tst == 0.) {
		goto L30;
	    }
	    if (tst <= sqrt((d__1 = D(m), abs(d__1))) * sqrt((d__2 = D(m + 1),
		     abs(d__2))) * eps) {
		E(m) = 0.;
		goto L30;
	    }
/* L20: */
	}
    }
    m = *n;

L30:
    l = l1;
    lsv = l;
    lend = m;
    lendsv = lend;
    l1 = m + 1;
    if (lend == l) {
	goto L10;
    }

/*     Scale submatrix in rows and columns L to LEND */

    i__1 = lend - l + 1;
    anorm = dlanst_("I", &i__1, &D(l), &E(l));
    iscale = 0;
    if (anorm == 0.) {
	goto L10;
    }
    if (anorm > ssfmax) {
	iscale = 1;
	i__1 = lend - l + 1;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &D(l), n, 
		info);
	i__1 = lend - l;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &E(l), n, 
		info);
    } else if (anorm < ssfmin) {
	iscale = 2;
	i__1 = lend - l + 1;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &D(l), n, 
		info);
	i__1 = lend - l;
	dlascl_("G", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &E(l), n, 
		info);
    }

/*     Choose between QL and QR iteration */

    if ((d__1 = D(lend), abs(d__1)) < (d__2 = D(l), abs(d__2))) {
	lend = lsv;
	l = lendsv;
    }

    if (lend > l) {

/*        QL Iteration   

          Look for small subdiagonal element. */

L40:
	if (l != lend) {
	    lendm1 = lend - 1;
	    i__1 = lendm1;
	    for (m = l; m <= lendm1; ++m) {
/* Computing 2nd power */
		d__2 = (d__1 = E(m), abs(d__1));
		tst = d__2 * d__2;
		if (tst <= eps2 * (d__1 = D(m), abs(d__1)) * (d__2 = D(m + 1),
			 abs(d__2)) + safmin) {
		    goto L60;
		}
/* L50: */
	    }
	}

	m = lend;

L60:
	if (m < lend) {
	    E(m) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L80;
	}

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2   
          to compute its eigensystem. */

	if (m == l + 1) {
	    if (icompz > 0) {
		dlaev2_(&D(l), &E(l), &D(l + 1), &rt1, &rt2, &c, &s);
		WORK(l) = c;
		WORK(*n - 1 + l) = s;
		dlasr_("R", "V", "B", n, &c__2, &WORK(l), &WORK(*n - 1 + l), &
			Z(1,l), ldz);
	    } else {
		dlae2_(&D(l), &E(l), &D(l + 1), &rt1, &rt2);
	    }
	    D(l) = rt1;
	    D(l + 1) = rt2;
	    E(l) = 0.;
	    l += 2;
	    if (l <= lend) {
		goto L40;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

/*        Form shift. */

	g = (D(l + 1) - p) / (E(l) * 2.);
	r = dlapy2_(&g, &c_b10);
	g = D(m) - p + E(l) / (g + d_sign(&r, &g));

	s = 1.;
	c = 1.;
	p = 0.;

/*        Inner loop */

	mm1 = m - 1;
	i__1 = l;
	for (i = mm1; i >= l; --i) {
	    f = s * E(i);
	    b = c * E(i);
	    dlartg_(&g, &f, &c, &s, &r);
	    if (i != m - 1) {
		E(i + 1) = r;
	    }
	    g = D(i + 1) - p;
	    r = (D(i) - g) * s + c * 2. * b;
	    p = s * r;
	    D(i + 1) = g + p;
	    g = c * r - b;

/*           If eigenvectors are desired, then save rotations. */

	    if (icompz > 0) {
		WORK(i) = c;
		WORK(*n - 1 + i) = -s;
	    }

/* L70: */
	}

/*        If eigenvectors are desired, then apply saved rotations. */

	if (icompz > 0) {
	    mm = m - l + 1;
	    dlasr_("R", "V", "B", n, &mm, &WORK(l), &WORK(*n - 1 + l), &Z(1,l), ldz);
	}

	D(l) -= p;
	E(l) = g;
	goto L40;

/*        Eigenvalue found. */

L80:
	D(l) = p;

	++l;
	if (l <= lend) {
	    goto L40;
	}
	goto L140;

    } else {

/*        QR Iteration   

          Look for small superdiagonal element. */

L90:
	if (l != lend) {
	    lendp1 = lend + 1;
	    i__1 = lendp1;
	    for (m = l; m >= lendp1; --m) {
/* Computing 2nd power */
		d__2 = (d__1 = E(m - 1), abs(d__1));
		tst = d__2 * d__2;
		if (tst <= eps2 * (d__1 = D(m), abs(d__1)) * (d__2 = D(m - 1),
			 abs(d__2)) + safmin) {
		    goto L110;
		}
/* L100: */
	    }
	}

	m = lend;

L110:
	if (m > lend) {
	    E(m - 1) = 0.;
	}
	p = D(l);
	if (m == l) {
	    goto L130;
	}

/*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2   
          to compute its eigensystem. */

	if (m == l - 1) {
	    if (icompz > 0) {
		dlaev2_(&D(l - 1), &E(l - 1), &D(l), &rt1, &rt2, &c, &s);
		WORK(m) = c;
		WORK(*n - 1 + m) = s;
		dlasr_("R", "V", "F", n, &c__2, &WORK(m), &WORK(*n - 1 + m), &
			Z(1,l-1), ldz);
	    } else {
		dlae2_(&D(l - 1), &E(l - 1), &D(l), &rt1, &rt2);
	    }
	    D(l - 1) = rt1;
	    D(l) = rt2;
	    E(l - 1) = 0.;
	    l += -2;
	    if (l >= lend) {
		goto L90;
	    }
	    goto L140;
	}

	if (jtot == nmaxit) {
	    goto L140;
	}
	++jtot;

/*        Form shift. */

	g = (D(l - 1) - p) / (E(l - 1) * 2.);
	r = dlapy2_(&g, &c_b10);
	g = D(m) - p + E(l - 1) / (g + d_sign(&r, &g));

	s = 1.;
	c = 1.;
	p = 0.;

/*        Inner loop */

	lm1 = l - 1;
	i__1 = lm1;
	for (i = m; i <= lm1; ++i) {
	    f = s * E(i);
	    b = c * E(i);
	    dlartg_(&g, &f, &c, &s, &r);
	    if (i != m) {
		E(i - 1) = r;
	    }
	    g = D(i) - p;
	    r = (D(i + 1) - g) * s + c * 2. * b;
	    p = s * r;
	    D(i) = g + p;
	    g = c * r - b;

/*           If eigenvectors are desired, then save rotations. */

	    if (icompz > 0) {
		WORK(i) = c;
		WORK(*n - 1 + i) = s;
	    }

/* L120: */
	}

/*        If eigenvectors are desired, then apply saved rotations. */

	if (icompz > 0) {
	    mm = l - m + 1;
	    dlasr_("R", "V", "F", n, &mm, &WORK(m), &WORK(*n - 1 + m), &Z(1,m), ldz);
	}

	D(l) -= p;
	E(lm1) = g;
	goto L90;

/*        Eigenvalue found. */

L130:
	D(l) = p;

	--l;
	if (l >= lend) {
	    goto L90;
	}
	goto L140;

    }

/*     Undo scaling if necessary */

L140:
    if (iscale == 1) {
	i__1 = lendsv - lsv + 1;
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &D(lsv), n, 
		info);
	i__1 = lendsv - lsv;
	dlascl_("G", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &E(lsv), n, 
		info);
    } else if (iscale == 2) {
	i__1 = lendsv - lsv + 1;
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &D(lsv), n, 
		info);
	i__1 = lendsv - lsv;
	dlascl_("G", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &E(lsv), n, 
		info);
    }

/*     Check for no convergence to an eigenvalue after a total   
       of N*MAXIT iterations. */

    if (jtot < nmaxit) {
	goto L10;
    }
    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
	if (E(i) != 0.) {
	    ++(*info);
	}
/* L150: */
    }
    goto L190;

/*     Order eigenvalues and eigenvectors. */

L160:
    if (icompz == 0) {

/*        Use Quick Sort */

	dlasrt_("I", n, &D(1), info);

    } else {

/*        Use Selection Sort to minimize swaps of eigenvectors */

	i__1 = *n;
	for (ii = 2; ii <= *n; ++ii) {
	    i = ii - 1;
	    k = i;
	    p = D(i);
	    i__2 = *n;
	    for (j = ii; j <= *n; ++j) {
		if (D(j) < p) {
		    k = j;
		    p = D(j);
		}
/* L170: */
	    }
	    if (k != i) {
		D(k) = D(i);
		D(i) = p;
		dswap_(n, &Z(1,i), &c__1, &Z(1,k), &
			c__1);
	    }
/* L180: */
	}
    }

L190:
    return 0;

/*     End of DSTEQR */

} /* dsteqr_ */

