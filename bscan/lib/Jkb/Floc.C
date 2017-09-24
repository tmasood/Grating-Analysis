#include <fstream>
#include <iostream>
#include <cstdlib>
#include <math.h>

using namespace std;

/*****************************************************************
  subroutine zloc

  purpose
     locating the zeros of a function of
     one variable in an interval

  usage
     call zloc(array,sub1,sub2)

  description
     array - array of length 5
            array(0) - endpoint of search interval
            array(1) - other endpoint
            array(2) - tentative step size
            array(3) - zero counter
            array(4) - termination switch
     sub1 - user supplied subroutine subr1(x,y)
            zloc supplies argument x and user stores
            corresponding function value in y
     sub2 - user supplied subroutine subr2(a,z,e)
            current zero z and associated absolute
            error e is supplied to user

  remarks
     1. array(3) contains as a floating number the
           current number of zeros found
     2. user may end run at any time prior to normal
           termination by making array(4) different from 0.
     3. search proceeds from array(0) to array(1)
        regardless fo the ordering of these entries
     4. an error stop with an appropriate message occurs is
        a - the search interval becomes too small
        b - 30 successive interval halvings occur
        c - the endpoints are initially equal

  method
    mueller,s method coupled with adaptive search
    based on quadratic extrapolation
*******************************************************************/
void Floc(double *array, double(*apzoo1)(double ), 
          void(*apzoo2)(double *, double, double ) )
{
  const int TRUE = 1;
  int Loop;
  int i, im, inst, ip, j, j1, j2, krt, nc, nd;
  double a, ah, alpha, b, beta, bh, ch, del, dh, dlamb, e; 
  double h, hpr, q0, q1, q2, r, rho, w, xp, yg, yp, z;
  double f[3], s[4], t[4], u[3], x[4], y[4];
  double loge;
 

  e = 1.e-15;
  loge = log( 1.e-15 )*0.05 - 0.0346,2;
  nd = 0;
  nc = 0;
  *(array+3) = 0;
  *(array+4) = 0;
  a = *array;
  b = *(array+1);
  h = *(array+2);
  if( h + a == a ){
    h = (b - a)/4;
  }
  else{
    h = fabs( h );
    if( h > fabs( b - a )/4)
      h = fabs( b - a )/4;
    if( b < a )
      h = -h;
  }
  if( h == 0 ){
    cerr <<  "Subroutine floc called with = endpoints" << endl;
    exit( EXIT_FAILURE );
  }
  else{
    z = a;
    alpha = a;
    beta = b;
    if( b < a )
      alpha = b;
    if( a > b )
      beta = a;
    j = 1;
    while( j < 5 ){
      w = (*apzoo1)( z );
      if( *(array+4) != 0 )
	goto L_250;
      *(x+j-1) = z;
      *(y+j-1) = w;
      z += h;
      if( alpha > z )
	z = alpha;
      if( beta < z )
	z = beta;
      if( z == *(x+j-1) )
	goto L_260;
      h = z - *(x+j-1);
      j++;
    }
    Loop = 1;
    while( Loop ){
      yg = *y + 3*( *(y+2) - *(y+1) );
      dh = fabs( *(y+3) );
      if( fabs( *(y+2) ) > dh )
	dh = fabs( *(y+2) );
      if( fabs( *(y+1) ) > dh )
	dh = fabs( *(y+1) );
      if( fabs( *y ) > dh )
	dh = fabs( *y );
      dh *= 3;
      if( fabs( *(y+3) - yg ) - 0.0333*dh > 0 ){
	nc += 1;
	if( nc >= 30 ){
	  cerr << " Search interval halved 30 times.\n" << endl;
	  exit ( EXIT_FAILURE );
	}
	h /= 2;
	*(x+2) = *(x+1);
	*(y+2) = *(y+1);
	z = a + h;
	j = 2;
	while( j < 5 ){
	  w = (*apzoo1)( z );
	  if( *(array+4) != 0 )
	    goto L_250;
	  *(x+j-1) = z;
	  *(y+j-1) = w;
	  z += h + h;
	  j += 2;
	}
	if( *x == *(x+1) || *(x+1) == *(x+2) || *(x+2) == *(x+3) ){
	  cerr << " Search interval is too small.\n" << endl;
	  exit ( EXIT_FAILURE );
	}
      }
      else
	Loop = 0;
    }
    nc = 0;
    inst = 1;
    i = 1;
    z = *x;
    w = *y;
    krt = -1;
    if( w == 0 )
      goto L_150;
  }
 L_100:
  krt = 0;
  z = *(x+i);
  w = *(y+i);
  if( w == 0 )
    goto L_150;
 L_110:
  krt = 1;
  if( *(y+i-1) == 0 )
    goto L_210;
  if( *(y+i-1)/fabs( *(y+i-1) )**(y+i) >= 0 )
    goto L_210;
  dlamb = x[i-1];
  del = y[i-1];
  rho = x[i];
  r = y[i];
  s[0] = x[0];
  s[1] = x[1];
  s[2] = x[2];
  t[0] = y[0];
  t[1] = y[1];
  t[2] = y[2];
 L_120:
  q2 = t[2]*(s[1] - s[0]) + t[0]*(s[2] - s[1]) + t[1]*(s[0] - s[2]);
  q1 = (t[1] - t[2])*s[0]*s[0];
  q1 += (t[2] - t[0])*s[1]*s[1];
  q1 = ((t[0] - t[1])*s[2]*s[2] + q1)/2;
  q0 = (s[1] - s[0])*s[0]*s[1]*t[2];
  q0 += (s[2] - s[1])*s[1]*s[2]*t[0];
  q0 += (s[0] - s[2])*s[2]*s[0]*t[1];
  if( q2 != 0 ){
    if( q1 == 0 ){
      z = -(q0/q2);
      if( z < 0 )
	goto L_130;
      if( z != 0 )
	z = sqrt(z);
    }
    else{
      z = fabs(q1) - (q0/fabs(q1))*q2;
      if( z < 0 )
	goto L_130;
      z = (fabs(q1) + sqrt(fabs(q1))*sqrt(z))/q2;
      if( q1 < 0 )
	z = -z;
    }
    if( (z - dlamb)*(rho - z) >= 0 )
      goto L_140;
    if( z == 0.e0 )
      goto L_130;
    z = q0/(q2*z);
  }
  else if( q1 == 0.e0 ){
    goto L_130;
  }
  else{
    z = q0/(q1*2);
  }
  if( (z - dlamb)*(rho - z) > 0.e0 )
    goto L_140;
 L_130:
  z = (r*dlamb - del*rho)/(r - del);
  if( (z - dlamb)*(rho - z) < 0.e0 )
    z = s[2];
 L_140:
  s[3] = z;
  w = (*apzoo1)( z );
  if( array[4] != 0.e0 )
    goto L_250;
  t[3] = w;
  if( w != 0.e0 ){
    if( fabs( z - s[2] ) - e*fabs( s[0] ) > 0.e0 ){
      q0 = fabs( w );
      if( del < 0.e0 )
	q0 = -q0;
      if( q0 == w ){
	dlamb = z;
	del = w;
      }
      else{
	rho = z;
	r = w;
      }
      s[0] = s[1];
      s[1] = s[2];
      s[2] = s[3];
      t[0] = t[1];
      t[1] = t[2];
      t[2] = t[3];
      goto L_120;
    }
  }
 L_150:
  s[0] = x[0];
  t[0] = y[0];
  s[1] = x[2];
  t[1] = y[2];
  if( z != x[1] ){
    if( z == x[0] ){
      s[0] = x[1];
      t[0] = y[1];
    }
    else{
      s[1] = x[1];
      t[1] = y[1];
    }
  }
  w = (s[1] - s[0])/2;
  s[0] = (s[0] - z)/w;
  s[1] = (s[1] - z)/w;
  t[2] = (fabs( t[0] ) + fabs( t[1] ))/w;
  if( t[2] == 0 )
    w *= ((fabs( s[0] ) + fabs( s[1] ))/2 + 1)/2;
  if( t[2] != 0 ){
    t[0] = t[0]/s[0]/t[2];
    t[1] = t[1]/s[1]/t[2];
    s[2] = fabs( t[1] - t[0] )*2/e;
    t[2] = fabs( t[1]*s[0] - t[0]*s[1] )/e;
    if( s[2] + t[2] != 0 ){
      s[2] = w/(t[2] + sqrt( fabs( t[2]*t[2] - s[2] ) ))*loge*loge;
      if( fabs( s[2] ) < fabs( w ) )
	w = s[2];
    }
  }
  if( fabs( w ) < fabs( z )*e/2 )
    w = z*e/2;
  array[3] += 1;
  (*apzoo2)( array, z, fabs(w) );
  if( array[4] != 0 )
    goto L_250;
  if( nd == 0 )
    goto L_200;
  nd -= 1;
  if( nd == 0 )
    goto L_230;
  if( f[1]/(*f + f[2]) <= 0 )
    goto L_170;
  u[1] = z;
  f[1] = 0;
 L_160:
  nd = 2;
  if( f[2] == 0 ){
    nd = 1;
    if( u[2] != b )
      goto L_230;
  }
  else if( *f == 0 ){
    nd = 1;
    if( u[0] != a )
      goto L_230;
  }
  else{
    if( fabs( *f ) > fabs( f[1] ) ){
      if( fabs( f[2] ) <= fabs( f[1] ) ){
	if( u[2] != b )
	  goto L_230;
      }
    }
    else if( u[0] != a ){
      goto L_230;
    }
    if( f[1] == 0 )
      nd = 1;
    if( *f/f[2] < 0 )
      goto L_230;
  }
  hpr = u[1];
  ch = f[1];
  while( TRUE ){
    if( fabs( f[1] ) < fabs( ch ) )
      hpr = u[1];
    if( fabs( f[1] ) < fabs( ch ) )
      ch = f[1];
    if( (hpr - u[0])*(u[2] - hpr) < 0 )
      goto L_230;
    while( TRUE ){
      im = 1;
      ip = 3;
      if( fabs( f[2] ) > fabs( *f ) ){
	ip = 1;
	im = 3;
      }
      if( fabs( f[1] ) >= fabs( f[ip - 1] ) )
	break;
      z = u[ip-1] - u[1];
      if( fabs( z ) >= fabs( u[1] - u[im-1] ) )
	break;
      z += u[im-1];
      w = (*apzoo1)( z );
      if( array[4] != 0 )
	goto L_250;
      if( (w - f[1])/(*f + f[2]) <= 0 ){
	f[ip - 1] = f[1];
	u[ip - 1] = u[1];
	u[1] = z;
	f[1] = w;
      }
      else if( (w - f[ip - 1])/(*f + f[2]) < 0 ){
	u[1] = z;
	f[1] = w;
      }
      else{
	u[im - 1] = z;
	f[im - 1] = w;
      }
    }
    t[0] = (*f - f[1])/(u[0] - u[1]);
    t[2] = (f[1] - f[2])/(u[1] - u[2]);
    s[0] = t[0]/(u[2] - u[0]);
    s[2] = t[2]/(u[2] - u[0]);
    z = (t[2] - t[0])/(u[2] - u[0])*2;
    if( (*f + f[2])/fabs( *f + f[2] )*z <= 0 )
      goto L_230;
    z = (s[2]*(u[0] - u[1]) + s[0]*(u[1] - u[2]))/z + u[1];
    if( (z - u[0])*(u[2] - z) <= 0 )
      goto L_230;
    if( z == hpr )
      goto L_230;
    if( z == u[1] )
      goto L_230;
    yg = (s[2]*(z - u[0]) - s[0]*(z - u[2]))*(z - u[1]) + f[1];
    if( yg/(*f + f[2]) >= 0.05 )
      goto L_230;
    w = (*apzoo1)( z );
    if( array[4] != 0 )
      goto L_250;
    if( w == 0 )
      goto L_150;
    if( w/(*f + f[2]) >= 0.05 )
      goto L_230;
    if( w/(*f + f[2]) <= 0 )
      goto L_180;
    if( fabs( yg - w ) <= 0.05*(fabs( w ) + fabs( yg )) )
      goto L_230;
    im = 1;
    if( (z - u[0])*(u[1] - z) < 0 )
      im = 3;
    if( fabs( w ) > fabs( f[im - 1] ) )
      goto L_230;
    if( fabs( w ) <= fabs( f[1] ) ){
      ip = 4 - im;
      im = 2;
      f[ip - 1] = f[1];
      u[ip-1] = u[1];
    }
    f[im-1] = w;
    u[im-1] = z;
  }
 L_170:
  s[0] = z;
  t[0] = 0;
  s[1] = u[1];
  t[1] = f[1];
  s[2] = u[2];
  t[2] = f[2];
  z = u[1];
  w = f[1];
  goto L_190;
 L_180:
  t[0] = *f;
  t[1] = f[1];
  t[2] = f[2];
  s[0] = u[0];
  s[1] = u[1];
  s[2] = u[2];
  if( (z - s[0])*(s[1] - z) < 0 ){
    s[0] = s[1];
    t[0] = t[1];
  }
  else{
    s[2] = s[1];
    t[2] = t[1];
  }
  t[1] = w;
  s[1] = z;
  f[1] = w;
  u[1] = z;
  if( t[0] != 0 ){
    dlamb = s[0];
    del = t[0];
    rho = z;
    r = w;
    goto L_120;
  }
 L_190:
  dlamb = z;
  del = w;
  rho = s[2];
  r = t[2];
  goto L_120;
 L_200:
  if( krt < 0 )
    goto L_100;
  if( krt == 0 )
    goto L_110;
  goto L_220;
 L_210:
  if( i != 2 ){
    i = 2;
    goto L_100;
  }
  else if( fabs( y[0] ) + fabs( y[1] ) == 0 ){
    goto L_230;
  }
  else if( fabs( y[1] ) + fabs( y[2] ) == 0 ){
    goto L_230;
  }
  else if( y[0]/(y[1] + y[2]) < 0 ){
    goto L_230;
  }
  else if( fabs( y[2] ) + fabs( y[0] ) == 0 ){
    goto L_230;
  }
  else{
    u[0] = x[0];
    u[1] = x[1];
    u[2] = x[2];
    *f = y[0];
    f[1] = y[1];
    f[2] = y[2];
    goto L_160;
  }
 L_220:
  if( i != 2 ){
    i = 2;
    goto L_100;
  }
 L_230:
  nd = 0;
  if( inst != 1 ){
    if( x[2] == b )
      goto L_250;
    w = -yp;
    x[3] = x[2];
    z = 0;
    j = 0;
    while( TRUE ){
      j1 = j + 1 - ((j + 1)/3)*3 + 1;
      j2 = j + 2 - ((j + 2)/3)*3 + 1;
      z += fabs( y[j] );
      w += y[j]*(xp - x[j1-1])/(x[j]-x[j1-1])*(xp-x[j2-1])/
	(x[j]-x[j2-1]);
      j += 1;
      if( j > 2 )
	break;
    }
    h = (x[2] - x[1])*3/2;
    if( w != 0 ){
      ah = (x[2] - x[1]) + (x[2] - x[0]);
      bh = (x[2] - x[1])*(x[2] - x[0]);
      ch = 0.021517*(x[0] - xp)*(x[1] - xp)*(x[2] - xp)*z/
	fabs( w );
      ch *= 0.666666666;
      while( TRUE ){
	hpr = ((2*h + ah)*h*h + ch)/((3*h + 2*ah)*
				     h + bh);
	dh = fabs( hpr ) - fabs( h );
	if( dh >= 0 )
	  break;
	h = hpr;
	if( dh + 0.1*fabs( h ) >= 0 )
	  break;
      }
    }
    while( TRUE ){
      z = x[2] + h;
      if( alpha > z )
	z = alpha;
      if( beta < z )
	z = beta;
      if( z == x[2] )
	goto L_260;
      if( z == x[3] )
	goto L_260;
      h = z - x[2];
      w = (*apzoo1)( z );
      if( array[4] != 0 )
	goto L_250;
      x[3] = z;
      y[3] = w;
      j = 0;
      while( j < 3 ){
	j1 = j + 1 - ((j + 1)/3)*3 + 1;
	j2 = j + 2 - ((j + 2)/3)*3 + 1;
	w += -y[j]*(z - x[j1-1])/(x[j]-x[j1-1])*(z - 
						 x[j2-1])/(x[j]-x[j2-1]);
	j += 1;
      }
      dh = fabs( y[3] );
      if( fabs( y[2] ) > dh )
	dh = fabs( y[2] );
      if( fabs( y[1] ) > dh )
	dh = fabs( y[1] );
      if( fabs( y[0] ) > dh )
	dh = fabs( y[0] );
      dh *= 3;
      if( fabs( w ) - 0.0333*dh <= 0 )
	break;
      nc += 1;
      if( nc >= 30 )
	goto L_240;
      h = (z - x[2])/2;
    }
    nc = 0;
  }
  inst = 0;
  xp = x[0];
  yp = y[0];
  x[0] = x[1];
  x[1] = x[2];
  x[2] = x[3];
  y[0] = y[1];
  y[1] = y[2];
  y[2] = y[3];
  goto L_100;
 L_240:
  cerr << "30 halves in subr zloc\n";
  exit( EXIT_FAILURE );
 L_250:
  if( h > 0 ){
    a = alpha;
    b = beta;
  }
  else{
    b = alpha;
    a = beta;
  }
  return;
 L_260:
  cerr << " interval too small in subr zloc\n" ;
  exit( EXIT_FAILURE );
} // end of function 
