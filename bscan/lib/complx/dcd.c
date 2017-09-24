#include <stdio.h>
#include <stdlib.h>
#include <complx.h>

complx dcd( double r, complx l )   /* retn the quotent of double and  complx */
{
   complx dz;
   double dr, di, r2;

   dr = l.re;
   di = l.im;
   r2 = dr*dr+di*di;

   if( r2 == 0.){
     fprintf(stderr, "Cannot divide by 0.0 \n");
     exit( EXIT_FAILURE );
   }

   dz.re =  dr*r/r2;
   dz.im = -di*r/r2;
   return( dz );
}
