#include <stdio.h>
#include <math.h>
#include <complx.h>
#include <float.h>

/****************************************************************
 * casin(z) calculates the complx trigonometric arc sine of z.
 * The result is in units of radians, and the real part is in the
 * first or fourth quadrant.
 * Author: Fullerton, W., (LANL)
 ***************************************************************/
complx casin( complx zinp )
{
  static int first = 1;
  static complx CI = { 0., 1. };

  static int nterms;
  static double rmin;

  int i;
  double r, twoi;
  complx z, z2, sqzp1, casin_v;

  if( first ){
    nterms = -0.4343*log(DBL_EPSILON);
    rmin = sqrt( 6.0*DBL_EPSILON );
    first = 0;
  }

  z = zinp;
  r = zabs(z);

  if( r > 0.1 ){
    if( zinp.re < 0.0 )
      z = cneg( zinp );
    sqzp1 = csqrt( cfadd( z, 1.0 ) );
    if( sqzp1.im < 0.)
      sqzp1 = cneg( sqzp1 );
    casin_v = fcsub(M_PI_2,cmul( CI, zlog(cadd(z,cmul(sqzp1,csqrt(cfsub( z, 1.0 )))))));
    if( casin_v.re > M_PI_2 )
      casin_v = fcsub( M_PI, casin_v );
    if( casin_v.re <= -M_PI_2  )
      casin_v = fcsub( -M_PI, casin_v );
    if( zinp.re < 0. )
      casin_v = cneg( casin_v );
  }
  else{
    casin_v = z;
    if( r >= rmin ){
      casin_v.re = 0.;
      casin_v.im = 0.;
      z2 = cmul(z,z);
      for( i = 1; i <= nterms; i++){
	twoi = (double )(2*(nterms - i) + 1);
	casin_v = fcadd(1.0/twoi, fcm(twoi,cdd(cmul(casin_v,z2),twoi+1.0)));
      }
      casin_v = cmul(z,casin_v);
    }
  }
  return( casin_v );
}

/**********************************************************/
complx cacos( complx z )
{
  complx w;

  w = casin( z );
  w.re = M_PI_2  -  w.re;
  w.im = -w.im;
  return( w );
}

















