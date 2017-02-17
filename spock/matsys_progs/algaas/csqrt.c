#include <stdio.h>
#include "util.h"

double fcabs(double real, double imag)
{
  double temp;
 
  if(real < 0)
    {
      real = -real;
    }
  if(imag < 0)
    {
      imag = -imag;
    }
  if(imag > real)
    {
      temp = real;
      real = imag;
      imag = temp;
    }
  if((real+imag) == real)
    {
        return(real);
    }
 
  temp = imag/real;
  temp = real*sqrt(1.0 + temp*temp);  /*overflow!!*/
  return(temp);
}


void csqrt(complex *r, complex *z)
{
  double mag, t;
  double zi = z->i, zr = z->r;
 
  if( (mag = fcabs(zr, zi)) == 0.)
    {
      r->r = r->i = 0.;
    }
  else if(zr > 0)
    {
      r->r = t = sqrt(0.5 * (mag + zr) );
      t = zi / t;
      r->i = 0.5 * t;
    }
  else
    {
      t = sqrt(0.5 * (mag - zr) );
      if(zi < 0)
	t = -t;
      r->i = t;
      t = zi / t;
      r->r = 0.5 * t;
    }
}
