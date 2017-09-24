#include <complx.h>

complx zconj(complx a)
{
   complx v;

   v.re = a.re;
   v.im = -a.im;

   return( v );
}
