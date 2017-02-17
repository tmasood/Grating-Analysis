/*
 * local_complex.c
 * a few complex functions
 *
 */
#include <math.h>
#include "local_complex.h"

/*
 *  convert complex to polar
 */
void z2polar(dcomplex z, double *r, double *phi) {

  *r = sqrt(__real__ z * __real__ z + __imag__ z * __imag__ z) ;
  if ( *r < 1e-50  ) {
    *phi = 0.0;
  }
  else {
    *phi = asin( __imag__ z / (*r) );
    if ( __real__ z < 0 ) {
      if ( __imag z > 0 ) {
        *phi = M_PI - *phi;
      }
      else {
        *phi = -M_PI - *phi;
      }
    }
  }
} /* z2polar */

/*
 *  convert polar to complex
 */
void polar2z(dcomplex *z, double r, double phi) {
  __real__ *z = r*cos(phi);
  __imag__ *z = r*sin(phi);
} /* z2polar */

/*
 *   abs( z )
 */
double zabs( dcomplex z ) {
  double z1;
  z1 =  __real__ z * __real__ z + __imag__ z * __imag__ z;
  return sqrt(z1);
} /* zabs */


double zabsq( dcomplex z )
{
  double z1;
  z1 =  zabs(z) * zabs(z);
  return sqrt(z1);
} /* zabsq */


/*
 *   arg( z )
 */
double arg( dcomplex z ) {
  double r, phi;
  z2polar(z, &r, &phi);
  return phi;
} /* arg */

/*
 *   exp( z )
 */
dcomplex zexp( dcomplex x ) {
  double z1;
  dcomplex z;
  z1 = exp( __real__ x );
  __real__ z = cos( __imag__ x )*z1;
  __imag__ z = sin( __imag__ x )*z1;
  return z;
} /* zexp */

/*
 *   exp( i x )
 */
dcomplex zexpR( double x ) {
  dcomplex z;
  __real__ z = cos( x );
  __imag__ z = sin( x );
  return z;
} /* zexp */

/*
 *   sin( z )
 */
dcomplex zsin( dcomplex x ) {
  dcomplex z1, z;
  /* multiply by i */
  __real__ z = - __imag__ x;
  __imag__ z = __real__ x;

  z1  = zexp(z);
  z1 -= zexp(-z);

  /* multiply by -i*0.5 */
  __real__ z = (__imag__ z1)*0.5;
  __imag__ z = - (__real__ z1)*0.5;

  return z;
} /* zsin */


/*
 *   cos( z )
 */
dcomplex zcos( dcomplex x ) {
  dcomplex z1, z;
  /* multiply by i */
  __real__ z = - __imag__ x;
  __imag__ z = __real__ x;

  z1  = zexp(z);
  z1 += zexp(-z);

  z = z1*0.5;

  return z;
} /* zsin */


/*
 *   sqrt( z )
 *   the branch is such that the result is in the 1st or 4th quadrant
 *   the root of a negative real is positive imaginary 
 */
dcomplex zsqrt( dcomplex z ) {
  double r, phi;
  dcomplex z2;

  if ( __imag__ z == 0.0 ) {
    if ( __real__ z >= 0.0 ) {
      __real__ z2 = sqrt(__real__ z);
      __imag__ z2 = 0.0;
    }
    else {
      __real__ z2 = 0.0;
      __imag__ z2 = sqrt(- __real__ z);
    }
    return z2;
  }
  z2polar(z, &r, &phi);
  r = sqrt(r);
  phi = 0.5*phi;
  polar2z(&z2, r, phi);
  return z2;
} /* zsqrt */ 


dcomplex zcnjg( dcomplex z )
{
  dcomplex r;
  double zi = __imag__ z;
  __real__ r = __real__ z;
  __imag__ r = -zi;
  return r;
}

