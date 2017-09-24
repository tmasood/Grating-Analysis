#include <iostream>
#include <cstdlib>

using namespace std;

/****************************************************************
 * Calculates the complex the ratio of Chebyshev Polynomials of 
 * the second
 * kind using the recurrence relation:
 *   rn(z) = 1/( 2*z - rn-1(z) )
 * where r0(z) = 1/2z
 *  * 
 * Author: Butler, Jerome (SMU)
 * August 1999
 ***************************************************************/

double dunounp1(int n, double z )
{
  int i;
  double x;
  double v;

  if( z == 0.0 ){
    cerr <<  "Division by zero in dunounp1.\n";
    exit( EXIT_FAILURE );
  }
  v = 0.0;
  for( i = 0 ; i < n+1 ; i++){
    x = 2.0*z - v;
    v = 1.0/x;
  }
  return( v );
}
