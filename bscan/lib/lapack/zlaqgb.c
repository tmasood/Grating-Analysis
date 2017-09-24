#include "f2c.h"

/* Subroutine */ int zlaqgb_(integer *m, integer *n, integer *kl, integer *ku,
	 doublecomplex *ab, integer *ldab, doublereal *r, doublereal *c, 
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    ZLAQGB equilibrates a general M by N band matrix A with KL   
    subdiagonals and KU superdiagonals using the row and scaling factors 
  
    in the vectors R and C.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    KL      (input) INTEGER   
            The number of subdiagonals within the band of A.  KL >= 0.   

    KU      (input) INTEGER   
            The number of superdiagonals within the band of A.  KU >= 0. 
  

    AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)   
            On entry, the matrix A in band storage, in rows 1 to KL+KU+1. 
  
            The j-th column of A is stored in the j-th column of the   
            array AB as follows:   
            AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)   

            On exit, the equilibrated matrix, in the same storage format 
  
            as A.  See EQUED for the form of the equilibrated matrix.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDA >= KL+KU+1.   

    R       (output) DOUBLE PRECISION array, dimension (M)   
            The row scale factors for A.   

    C       (output) DOUBLE PRECISION array, dimension (N)   
            The column scale factors for A.   

    ROWCND  (output) DOUBLE PRECISION   
            Ratio of the smallest R(i) to the largest R(i).   

    COLCND  (output) DOUBLE PRECISION   
            Ratio of the smallest C(i) to the largest C(i).   

    AMAX    (input) DOUBLE PRECISION   
            Absolute value of largest matrix entry.   

    EQUED   (output) CHARACTER*1   
            Specifies the form of equilibration that was done.   
            = 'N':  No equilibration   
            = 'R':  Row equilibration, i.e., A has been premultiplied by 
  
                    diag(R).   
            = 'C':  Column equilibration, i.e., A has been postmultiplied 
  
                    by diag(C).   
            = 'B':  Both row and column equilibration, i.e., A has been   
                    replaced by diag(R) * A * diag(C).   

    Internal Parameters   
    ===================   

    THRESH is a threshold value used to decide if row or column scaling   
    should be done based on the ratio of the row or column scaling   
    factors.  If ROWCND < THRESH, row scaling is done, and if   
    COLCND < THRESH, column scaling is done.   

    LARGE and SMALL are threshold values used to decide if row scaling   
    should be done based on the absolute size of the largest matrix   
    element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1;
    doublecomplex z__1;
    /* Local variables */
    static integer i, j;
    static doublereal large, small, cj;
    extern doublereal dlamch_(char *);


#define R(I) r[(I)-1]
#define C(I) c[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    if (*m <= 0 || *n <= 0) {
	*(unsigned char *)equed = 'N';
	return 0;
    }

/*     Initialize LARGE and SMALL. */

    small = dlamch_("Safe minimum") / dlamch_("Precision");
    large = 1. / small;

    if (*rowcnd >= .1 && *amax >= small && *amax <= large) {

/*        No row scaling */

	if (*colcnd >= .1) {

/*           No column scaling */

	    *(unsigned char *)equed = 'N';
	} else {

/*           Column scaling */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		cj = C(j);
/* Computing MAX */
		i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
		i__5 = *m, i__6 = j + *kl;
		i__4 = min(i__5,i__6);
		for (i = max(1,j-*ku); i <= min(*m,j+*kl); ++i) {
		    i__2 = *ku + 1 + i - j + j * ab_dim1;
		    i__3 = *ku + 1 + i - j + j * ab_dim1;
		    z__1.r = cj * AB(*ku+1+i-j,j).r, z__1.i = cj * AB(*ku+1+i-j,j).i;
		    AB(*ku+1+i-j,j).r = z__1.r, AB(*ku+1+i-j,j).i = z__1.i;
/* L10: */
		}
/* L20: */
	    }
	    *(unsigned char *)equed = 'C';
	}
    } else if (*colcnd >= .1) {

/*        Row scaling, no column scaling */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
	    i__5 = *m, i__6 = j + *kl;
	    i__3 = min(i__5,i__6);
	    for (i = max(1,j-*ku); i <= min(*m,j+*kl); ++i) {
		i__4 = *ku + 1 + i - j + j * ab_dim1;
		i__2 = i;
		i__5 = *ku + 1 + i - j + j * ab_dim1;
		z__1.r = R(i) * AB(*ku+1+i-j,j).r, z__1.i = R(i) * AB(*ku+1+i-j,j).i;
		AB(*ku+1+i-j,j).r = z__1.r, AB(*ku+1+i-j,j).i = z__1.i;
/* L30: */
	    }
/* L40: */
	}
	*(unsigned char *)equed = 'R';
    } else {

/*        Row and column scaling */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    cj = C(j);
/* Computing MAX */
	    i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
	    i__5 = *m, i__6 = j + *kl;
	    i__2 = min(i__5,i__6);
	    for (i = max(1,j-*ku); i <= min(*m,j+*kl); ++i) {
		i__3 = *ku + 1 + i - j + j * ab_dim1;
		d__1 = cj * R(i);
		i__4 = *ku + 1 + i - j + j * ab_dim1;
		z__1.r = d__1 * AB(*ku+1+i-j,j).r, z__1.i = d__1 * AB(*ku+1+i-j,j).i;
		AB(*ku+1+i-j,j).r = z__1.r, AB(*ku+1+i-j,j).i = z__1.i;
/* L50: */
	    }
/* L60: */
	}
	*(unsigned char *)equed = 'B';
    }

    return 0;

/*     End of ZLAQGB */

} /* zlaqgb_ */

