#include <stdio.h>
#include <math.h>
#include <complx.h>

/****************************************************************
 * Calculates the complx Chebyshev Polynomials of the second
 * kind using the recurrence relation:
 *   Un+1(z) = 2*z*Un(z) - Un-1(z)
 * where U0(z) = 1, U1(z) = 2*z.
 * See Galler, Comm. ACM, June 1960.
 * 
 * Author: Butler, Jerome (SMU)
 ***************************************************************/

void zcheb(int n, complx z, complx *un, complx *unm1)
{
  int i;
  complx a, b, c;

  a.re = 1.0; a.im = 0.0;
  b.re = 2.0*z.re;
  b.im = 2.0*z.im;

  if( n == 0 ){
    c = a;
    b = a;
  }
  else{
    if( n == 1 ){
      c = b;
      b = a;
    }
    else{
      for( i = 2; i < n ; i++){
	c = csub(fcm(2.0,cmul(z,b)),a);
	a = b;
	b = c;
      }
      c = csub(fcm(2.0,cmul(z,b)),a);
    }
  }
  *un = c;
  *unm1 = b;
}
















