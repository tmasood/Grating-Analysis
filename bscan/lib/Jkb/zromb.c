#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complx.h"
#include "jkb.h"


complex zromb(eval, a, b, eps)
complex (*eval)();
double a, b, eps;
{
  int i, j, k, n;
  double d, h, x;
  complex endpts, sum, sumt, sumu, t[50], tprev[50], u[50], 
   uprev[50];

  endpts = fcm(0.5,cadd((*eval)(a),(*eval)(b)));
  sumt.re = 0.0, sumt.im = 0;
  i = 0;
  n = 1;
  h = b - a;
  while( 1 ){
    t[0] = fcm(h,cadd(endpts,sumt));
    sumu.re = 0.0, sumu.im = 0.0;
    x = a-0.5*h;
    for( j = 0; j < n; j++ ){
      x += h;
      sumu = cadd(sumu,(*eval)( x ));
      }
    u[0] = fcm(h,sumu);
    for(k = 0; k < i ; k++){
      if( cabs1( csub(t[k],u[k]) ) < eps ){
        sum = fcm(0.5,cadd(t[k],u[k]));
        return( sum );
        }
      else{
        d = 4.0*powi(4.0,k);
        t[k+1] = cfdiv(csub(fcm(d,t[k]),tprev[k]),(d-1.0));
        tprev[k] = t[k];
        u[k+1] = cfdiv(csub(fcm(d,u[k]),uprev[k]),(d-1.0));
        uprev[k] = u[k];
        }
      }
    h = 0.5*h;
    sumt = cadd(sumt,sumu);
    tprev[k] = t[k];
    uprev[k] = u[k];
    if( ++i >= 50 ){
      fprintf( stdout, "Error in function romb: faild to converge.\n" );
      exit(0);
      }
    n=n<<1;
    }
}

