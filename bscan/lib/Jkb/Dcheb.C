#include <stdio.h>
#include <math.h>

/****************************************************************
 * Calculates the complex Chebyshev Polynomials of the second
 * kind using the recurrence relation:
 *   Un+1(z) = 2*z*Un(z) - Un-1(z)
 * where U0(z) = 1, U1(z) = 2*z.
 * See Galler, Comm. ACM, June 1960.
 * 
 * Author: Butler, Jerome (SMU)
 ***************************************************************/

void Dcheb(int n, double z, double *un, double *unm1)
{
  int i;
  double a, b, c;

  a = 1.0;
  b = 2.0*z;

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
	c = 2*z*b-a;
	a = b;
	b = c;
      }
      c = 2*z*b-a;
    }
  }
  *un = c;
  *unm1 = b;
}
