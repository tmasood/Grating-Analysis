#include "f2c.h"

doublereal zlanht_(char *norm, integer *n, doublereal *d, doublecomplex *e)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1992   


    Purpose   
    =======   

    ZLANHT  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    complex Hermitian tridiagonal matrix A.   

    Description   
    ===========   

    ZLANHT returns the value   

       ZLANHT = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
  
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
  
    normF  denotes the  Frobenius norm of a matrix (square root of sum of 
  
    squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in ZLANHT as described   
            above.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, ZLANHT is   
            set to zero.   

    D       (input) DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of A.   

    E       (input) COMPLEX*16 array, dimension (N-1)   
            The (n-1) sub-diagonal or super-diagonal elements of A.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);
    /* Local variables */
    static integer i;
    static doublereal scale;
    extern logical lsame_(char *, char *);
    static doublereal anorm;
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *), zlassq_(integer *, doublecomplex *, 
	    integer *, doublereal *, doublereal *);
    static doublereal sum;



#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    if (*n <= 0) {
	anorm = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	anorm = (d__1 = D(*n), abs(d__1));
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
/* Computing MAX */
	    d__2 = anorm, d__3 = (d__1 = D(i), abs(d__1));
	    anorm = max(d__2,d__3);
/* Computing MAX */
	    d__1 = anorm, d__2 = z_abs(&E(i));
	    anorm = max(d__1,d__2);
/* L10: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1' || 
	    lsame_(norm, "I")) {

/*        Find norm1(A). */

	if (*n == 1) {
	    anorm = abs(D(1));
	} else {
/* Computing MAX */
	    d__2 = abs(D(1)) + z_abs(&E(1)), d__3 = z_abs(&E(*n - 1)) + (d__1 
		    = D(*n), abs(d__1));
	    anorm = max(d__2,d__3);
	    i__1 = *n - 1;
	    for (i = 2; i <= *n-1; ++i) {
/* Computing MAX */
		d__2 = anorm, d__3 = (d__1 = D(i), abs(d__1)) + z_abs(&E(i)) 
			+ z_abs(&E(i - 1));
		anorm = max(d__2,d__3);
/* L20: */
	    }
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	if (*n > 1) {
	    i__1 = *n - 1;
	    zlassq_(&i__1, &E(1), &c__1, &scale, &sum);
	    sum *= 2;
	}
	dlassq_(n, &D(1), &c__1, &scale, &sum);
	anorm = scale * sqrt(sum);
    }

    ret_val = anorm;
    return ret_val;

/*     End of ZLANHT */

} /* zlanht_ */

