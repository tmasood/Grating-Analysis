#include <stdio.h>
#include <math.h>
#include "complx.h"

complex  croot(complex (*fcn)(complex ), complex z1, complex z2, 
	       complex z3, complex *w4, int *i, int m, int *root)
{
   double r1, r2, r3;
   static double eps = 4.e-10;
   complex d1, d2, d3, om, om1, om2, w1, w2, w3 , z4;

   w1 = (*fcn)( z1 );
   w2 = (*fcn)( z2 );
   w3 = (*fcn)( z3 );
   *i = 3;
   d1 = cdiv((csub(w2,w1)),(csub(z2,z1)));
   while( *i < m ){
      d2 = cdiv((csub(w3,w2)),(csub(z3,z2)));
      d3 = cdiv((csub(d2,d1)),(csub(z3,z1)));
      om1 = cadd(d2,cmul((csub(z3,z2)),d3));
      om2 = zsqrt( csub(cmul(om1,om1),cmul(fcm(4.0,w3),d3)) );
      om = cadd(om1,om2);
      r1 = zabs( om );
      r2 = zabs( csub(om1,om2) );
      if( r1 < r2 )
         om = csub(om1,om2);
      z4 = csub(z3,cdiv(fcm(2,w3),om));
      *w4 = (*fcn)( z4 );
      r3 = zabs( csub(z4,z3) );
      if( r3 < eps ){
         *root = 1;
         return( z4 );
         }
      z1 = z2;
      w1 = w2;
      z2 = z3;
      w2 = w3;
      z3 = z4;
      w3 = *w4;
      d1 = d2;
      ++*i;
      }
   fprintf(stdout, "Failed to converge in %d iterations.\n", m );
   *root = 0;
   return( ftoc( 0.0, 0.0) );
   }
