#include <stdio.h>
#include <complx.h>

complx cdiv( complx l, complx r ) /* retn the quotient of two dp complx nos. */
{
   complx dz;
   double den;

   if( r.re == 0. && r.im == 0. ){
      fprintf( stderr, "complx division by 0.\n" );
      exit(0);
      }
   den = r.re*r.re + r.im*r.im;

   dz.re = ( l.re*r.re + l.im*r.im )/den;
   dz.im = ( r.re*l.im - l.re*r.im )/den;
   return( dz );
}
