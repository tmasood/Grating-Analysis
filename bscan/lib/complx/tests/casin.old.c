#include <stdio.h>
#include <math.h>
#include "complx.h"


complex casin( complex z )
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

/* Power series expansion */
/*
  b = cabs(z);
  if( b < 0.125 )
  {
  z2.re = (x - y) * (x + y);
  z2.im = 2.0 * x * y;

  cn = 1.0;
  n = 1.0;
  ca.re = x;
  ca.im = y;
  sum.re = x;
  sum.im = y;
  do
    {
    ct.re = z2.re * ca.re  -  z2.im * ca.im;
    ct.im = z2.re * ca.im  +  z2.im * ca.re;
    ca.re = ct.re;
    ca.im = ct.im;

    cn *= n;
    n += 1.0;
    cn /= n;
    n += 1.0;
    b = cn/n;
  
    ct.re *= b;
    ct.im *= b;
    sum.re += ct.re;
    sum.im += ct.im;
    b = fabs(ct.re) + fabs(ct.im);
    }
  while( b > MACHEP );
  w.re = sum.re;
  w.im = sum.im;
  return( w );
  }
  */
  

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
complex cacos( complex z )
{
  complex w;

  w = casin( z );
  w.re = 0.5*M_PI  -  w.re;
  w.im = -w.im;
  return( w );
}

