#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double powi(double, int );

double romb(eval, a, b, eps)
double (*eval)(), a, b, eps;
/************************************************************
 *  Function romb
 *
 *  Purpose
 *    Definite integration of a double function.
 *
 *  Parameters
 *    eval - double function to be integrated.
 *    a - lower limit.
 *    b - upper limit,
 *    eps - precision.
 *  
 *  Return Value
 *    Value of the integral.
 *
 *  Method
 *    Havie-Romberg Quadrature
 *    Reference: Algorithm 257 Comm ACM B(1965) 381
 *
 *  January 3, 1997
 *  Copyright c Jerome K. Butler
 ***********************************************************/
{
  int i, j, k, n;
  double d, endpts, h, sum, sumt, sumu, t[50], tprev[50], u[50], 
   uprev[50], x;

  endpts = 0.5*((*eval)(a)+(*eval)(b));
  sumt = 0.0;
  i = 0;
  n = 1;
  h = b - a;
  while( 1 ){
    t[0] = h*(endpts+sumt);
    sumu = 0.0;
    x = a-0.5*h;
    for( j = 0; j < n; j++ ){
      x += h;
      sumu += (*eval)( x );
      }
    u[0] = h*sumu;
    for(k = 0; k < i ; k++){
      if( fabs( t[k]-u[k] ) < eps ){
        sum = 0.5*(t[k]+u[k]);
        return( sum );
        }
      else{
        d = 4.0*powi(4.0,k);
        t[k+1] = (d*t[k]-tprev[k])/(d-1.0);
        tprev[k] = t[k];
        u[k+1] = (d*u[k]-uprev[k])/(d-1.0);
        uprev[k] = u[k];
        }
      }
    h = 0.5*h;
    sumt = sumt + sumu;
    tprev[k] = t[k];
    uprev[k] = u[k];
    if( ++i >= 50 ){
      fprintf( stdout, "Error in function romb: faild to converge.\n" );
      exit(0);
      }
    n=n<<1;
    }
}

