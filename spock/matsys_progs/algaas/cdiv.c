#include <stdio.h>
#include "util.h"

void cdiv(complex *c, complex *a, complex *b)
{
  double ratio, den;
  double abr, abi, cr;
  float af, bf;
  
  if( (abr = b->r) < 0.0)
    {
      abr = -abr;
    }

  if( (abi = b->i) < 0.0)
    {
      abi = -abi;
    }

  if( abr <= abi )
    {
      if(abi == 0)
	{
	  af = bf = abr;
	  if (a->i != 0 || a->r != 0)
	    {
	      af = 1.0;
	    }
	  c->i = c->r = af / bf;
	  return;
      }
      ratio = (double)b->r / b->i ;
      den = b->i * (1 + ratio*ratio);
      cr = (a->r*ratio + a->i) / den;
      c->i = (a->i*ratio - a->r) / den;
    }
  else
    {
      ratio = (double)b->i / b->r ;
      den = b->r * (1 + ratio*ratio);
      cr = (a->r + a->i*ratio) / den;
      c->i = (a->i - a->r*ratio) / den;
    }
  c->r = cr;
}
