#include <stdio.h>
#include <stdlib.h>
#include <complx.h>

complx cdd( complx r, double l )   /* retn the product of two dp complx nos. */
{
   complx dz;

   if( l == 0.){
     fprintf(stderr, "Cannot divide by 0.0 \n");
     exit( EXIT_FAILURE );
   }

   dz.re = r.re/l;
   dz.im = r.im/l;
   return( dz );
}
