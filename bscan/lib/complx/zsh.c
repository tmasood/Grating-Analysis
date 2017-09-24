#include <math.h>
#include <complx.h>

complx zsh( complx z)
{
   double ch, cn, sh, sn;
   complx zsh_v;

   sh = sinh(z.re);
   ch = cosh(z.re);
   sn = sin(z.im);
   cn = cos(z.im);
   zsh_v.re = sh*cn;
   zsh_v.im = ch*sn;
   return(zsh_v);
}
