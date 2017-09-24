#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complx.h>

/****************************************************************
 * Calculates the complx the ratio of Chebyshev Polynomials of 
 * the second
 * kind using the recurrence relation:
 *   rn(z) = 1/( 2*z - rn-1(z) )
 * where r0(z) = 1/2z
 *  * 
 * Author: Butler, Jerome (SMU)
 ***************************************************************/

complx zunounp1(int n, complx z )
{
  int i;
  double x, y, r2;
  complx v;

  if( (z.re == 0.0) && (z.im == 0) ){
    fprintf(stderr, "Division by zero in zunounp1.\n");
    exit( EXIT_FAILURE );
  }
  v.re = 0.0;
  v.im = 0.0;
  for( i = 0 ; i < n+1 ; i++){
    x = 2.0*z.re - v.re;
    y = 2.0*z.im - v.im;
    r2 = x*x + y*y;
    v.re = x/r2;
    v.im = -y/r2;
  }
  return( v );
}
















