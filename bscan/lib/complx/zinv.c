#include <complx.h>

complx zinv(complx a)
{
   double bm;
   complx v;

   bm = a.re*a.re+a.im*a.im;
   v.re = a.re/bm;
   v.im =  - a.im/bm;
   return( v );
}
