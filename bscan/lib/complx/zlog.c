#include <stdio.h>
#include <math.h>
#include <complx.h>

complx zlog( complx a )
{
   double dtheta, zm;
   complx zlog_v;

   if( a.re == 0.0 ){
      if( a.im == 0.0 ){
         printf( " zlog called with zero argument\n" );
         exit(0);
         }
      else{
         zlog_v.im = M_PI_2;
         zlog_v.re = log( fabs( a.im ) );
         if( a.im < 0.0 )
            zlog_v.im = -zlog_v.im;
         }
      }
   else if( a.im != 0.0 ){
      dtheta = atan( a.im/a.re );
      if( dtheta <= 0.0 ){
         if( a.re < 0.0 )
            dtheta += M_PI;
         }
      else{
         if( a.re < 0.0 )
            dtheta -= M_PI;
         }
      zm = zabs( a );
      zlog_v.re = log( zm );
      zlog_v.im = dtheta;
      }
   else if( a.re > 0.0 ){
      zlog_v.re = log( a.re );
      zlog_v.im = 0.0;
      }
   else{
      zlog_v.re = log( fabs( a.re ) );
      zlog_v.im = M_PI;
      }
   return( zlog_v );
}
