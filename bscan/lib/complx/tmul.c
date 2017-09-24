#include <complx.h>

complx tmul( complx l, complx c, complx r )   /* retn the product of three dp complx nos. */
{
   complx z;
   
   z.re = -r.im*(c.re*l.im+c.im*l.re)+r.re*(-c.im*l.im+c.re*l.re);
   z.im =  r.im*(-c.im*l.im+c.re*l.re)+r.re*(c.re*l.im+c.im*l.re);

   return( z );
}
