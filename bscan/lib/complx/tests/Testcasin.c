#include <stdio.h>
#include <math.h>
#include <complx.h>

complex casinold( complex z )
{
static complex ca, ct, zz, z2;
double x, y;
complex w;

  x = z.re;
  y = z.im;

  if( y == 0.0 )
    {
    if( fabs(x) > 1.0 )
      {
      w.re = 0.5*M_PI;
      w.im = 0.0;
      printf( "casin - DOMAIN ERROR\n" );
      }
    else
      {
      w.re = asin(x);
      w.im = 0.0;
      }
    return ( w );
    }

  ca.re = x;
  ca.im = y;

  ct.re = -ca.im;  /* iz */
  ct.im = ca.re;

  /* sqrt( 1 - z*z) */
/* cmul( &ca, &ca, &zz ) */

  zz.re = (ca.re - ca.im) * (ca.re + ca.im);  /*x * x  -  y * y */
  zz.im = 2.0 * ca.re * ca.im;

  zz.re = 1.0 - zz.re;
  zz.im = -zz.im;
  z2 = zsqrt( zz );

  zz = cadd( z2, ct );
  zz = zlog( zz );
  w.re = zz.im;  /* mult by 1/i = -i */
  w.im = -zz.re;
  return( w );
  }
/**********************************************************/
complex cacosold( complex z )
{
  complex w;

  w = casinold( z );
  w.re = 0.5*M_PI  -  w.re;
  w.im = -w.im;
  return( w );
}

/****************************************************************
 * casin(z) calculates the complex trigonometric arc sine of z.
 * The result is in units of radians, and the real part is in the
 * first or fourth quadrant.
 * Author: Fullerton, W., (LANL)
 ***************************************************************/
complex casinnew( complex zinp )
{
  static int first = 1;
  static complex CI = { 0., 1. };

  static int nterms;
  static double rmin;

  int i;
  double r, twoi;
  complex z, z2, sqzp1, casin_v;

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
	twoi = 2*(nterms - i) + 1;
	casin_v = fcadd(1.0/twoi, fcm(twoi,cfdiv(cmul(casin_v,z2),twoi+1.0)));
      }
      casin_v = cmul(z,casin_v);
    }
  }
  return( casin_v );
}

/**********************************************************/
complex cacosnew( complex z )
{
  complex w;

  w = casinnew( z );
  w.re = M_PI_2  -  w.re;
  w.im = -w.im;
  return( w );
}

int main()
{
  complex w, z;

  while( scanf("%lf %lf",&z.re,&z.im) > 1 ){
    printf(" z = %f  %f\n", z.re, z.im);
    w = casinnew( z );
    printf(" casinnew = %20.15f  %20.15f\n", w.re, w.im);
    w = cacosnew( z );
    printf(" cacosnew = %20.15f  %20.15f\n", w.re, w.im);

    w = casinold( z );
    printf(" casinold = %20.15f  %20.15f\n", w.re, w.im);
    w = cacosold( z );
    printf(" cacosold = %20.15f  %20.15f\n", w.re, w.im);
    }
  return( 0 );
}
