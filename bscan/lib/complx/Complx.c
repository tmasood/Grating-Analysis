#include <stdio.h>
#include <math.h>
#include <complx.h>

complx cneg(complx dz)   /* dp complx unary minus operation */
{
   dz.re = -dz.re;   dz.im = -dz.im;
   return( dz );
}

complx cadd(complx l, complx r)   /* retn the sum of two dp complx nos. */
{
   l.re += r.re;
   l.im += r.im;
   return( l );
}

complx csub( complx l, complx r )   /* retn the difference of two dp complx nos. */
{
   l.re -= r.re;
   l.im -= r.im;
   return( l );
}

complx cmul( complx l, complx r )   /* retn the product of two dp complx nos. */
{
   complx dz;
   dz.re = l.re*r.re - l.im*r.im;
   dz.im = l.re*r.im + l.im*r.re;
   return( dz );
}

complx tmul( complx l, complx c, complx r )   /* retn the product of three dp complx nos. */
{
   complx z;
   
   z = cmul(l,cmul(c,r)); 
   return( z );
}

complx fmul(a, b, c, d)
complx a, b, c, d;
{
   complx z;

   z = cmul(a,cmul(b,cmul(c,d)));
   return( z );
}

complx fcm( double l, complx r )   /* retn the product of two dp complx nos. */
{
   complx dz;
   dz.re = l*r.re;
   dz.im = l*r.im;
   return( dz );
}

complx cdiv( l, r )   /* retn the quotient of two dp complx nos. */
complx l, r;
{
   complx dz;
   double den;

   if( r.re == 0. && r.im == 0. ){
      printf( "complx division by 0." );
      exit(0);
      }
   den = r.re*r.re + r.im*r.im;

   dz.re = ( l.re*r.re + l.im*r.im )/den;
   dz.im = ( r.re*l.im - l.re*r.im )/den;
   return( dz );
}

complx tdiv( l, c, r )   /* retn l*c/r complx nos. */
complx l, c, r;
{
   complx dz;

   dz = cmul(l,cdiv(c,r));
   return( dz );
}

complx fdiv( a, b, c, d)
complx a, b, c, d;
{
   complx z;

   z = cdiv(tmul(a,b,c),d);
   return( z );
}   

complx dconjg( complx z)
{
   complx w;
   w.re = z.re;
   w.im =-z.im;
   return( w );
}

int ceq( l, r )   /* return TRUE if dp complx nos. equal */
complx l, r;
{
   if( l.re == r.re  &&  l.im == r.im )
      return( 1 );
   return( 0 );
}

complx ftoc( r, i )   /* convert doubles to dp complx */
double r, i;
{
   complx dz;
   dz.re = r;   dz.im = i;
   return( dz );
}

complx zsqrt( complx a )
{
   double dtheta, zm;
   complx v;

   zm = sqrt( sqrt( a.re*a.re+a.im*a.im ) );
   if( a.re == 0.0 ){           /**    The y axis */
      if( a.im > 0.0 ){         /**/
         v.re = zm*M_SQRT1_2;   /**    positive y axis */
         v.im = zm*M_SQRT1_2;   /**/
         }                      /**/
      else if( a.im < 0.0 ){    /**/
         v.re = zm*M_SQRT1_2;   /**    negative y axis */
         v.im = -zm*M_SQRT1_2;  /**/
         }                      /**/
      else{                     /**/
         v.re = 0.0;            /**    zero point */
         v.im = 0.0;            /**/
         }
      }
   else if( a.im != 0.0 ){        /* real part not zero */
      dtheta = atan( a.im/a.re );
      if( dtheta <= 0.0 ){        /***/ 
         if( a.re < 0.0 )         /***  second quadrant */
            dtheta += M_PI;       /***/
         }
      else{
         if( a.re < 0.0 )            /****  third quadrant */
            if( dtheta < 0.25*M_PI ) /****  br cut at 225 degrees */
               dtheta += M_PI;       /****/
            else                     /****/ 
               dtheta -= M_PI;       /****/
         }
      dtheta *= 0.5;                  /* first and forth quadrant */
      v.re = zm*cos( dtheta );  
      v.im = zm*sin( dtheta );
      }
   else if( a.re > 0.0 ){             /* positive x axis */
      v.re = sqrt( a.re );
      v.im = 0.0;
      }
   else{
      v.re = 0.0;               /* negative x axis */
      v.im = sqrt( fabs( a.re ) );
      }

   return( v );
}

complx zlog(a)
complx a;
{
   double dtheta, zabs(), zm;
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

complx zexp(a)
complx a;
{
   double zm;
   complx v;

   zm = exp( a.re );
   v.re = zm*cos(a.im);
   v.im = zm*sin(a.im);
   return( v );
}

complx csin( dz )   /* sine of dp complx no. */
complx dz;
{
   complx dzs;

   dzs.re = sin( dz.re )*cosh( dz.im );
   dzs.im = cos( dz.re )*sinh( dz.im );
   return( dzs );
}

complx ccos( dz )   /* cosine of dp complx no. */
complx dz;
{
   complx dzc;

   dzc.re = cos(dz.re)*cosh(dz.im);
   dzc.im = -sin(dz.re)*sinh(dz.im);
   return( dzc );
}

double cabs1(z)
complx z;
{
   double w;
   w = fabs(z.re) + fabs(z.im);
   return( w );
}   

complx ztan(complx z)
{
   double t1, u, v, x, y;
   complx w;

   x =  z.re;
   y =  z.im;
   u = tan(x);
   v = tanh(y);
   t1 = 1.0 + u*u*v*v;
   w.re = u*(1.0 - v*v)/t1;
   w.im = v*(1.0 + u*u)/t1;
   return ( w );
} 


double zabs(complx z)
{
   double q, s, u, v, zabs_v;

   u = fabs( z.re );
   v = fabs( z.im );
   s = u + v;
   s *= 1.0e0;
   if( s == 0.0e0 ){
      zabs_v = 0.0e0;
      }
   else if( u > v ){
      q = v/u;
      zabs_v = u*sqrt( 1.e0 + q*q );
      }
   else{
      q = u/v;
      zabs_v = v*sqrt( 1.e0 + q*q );
      }
   return( zabs_v );
}

complx zdiv(complx a, complx b)
{
   double bm, cc, cd, zabs();
   complx zdiv_v;

   bm = 1.0e0/zabs( b );
   cc = b.re*bm;
   cd = b.im*bm;
   zdiv_v.re = (a.re*cc + a.im*cd)*bm;
   zdiv_v.im = (a.im*cc - a.re*cd)*bm;
   return( zdiv_v );
}

complx zch(complx z)
{
   double ch, cn, sh, sn;
   complx zch_v;

   sh = sinh(z.re);
   ch = cosh(z.re);
   sn = sin(z.im);
   cn = cos(z.im);
   zch_v.re = ch*cn;
   zch_v.im = sh*sn;
   return( zch_v );
}

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

/*=================================================================== */
complx  croot(fcn, z1, z2, z3, w4, i, m, root)
complx (*fcn)(complx ), z1, z2, z3, *w4;
int *i, m, *root;
{
   double r1, r2, r3;
   static double eps = 4.e-10;
   complx d1, d2, d3, om, om1, om2, w1, w2, w3 , z4;

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
/*=================================================================== */
void cguesses(complx *z1, complx *z2, complx z3)
{
   *z1 = cadd(cmul(ftoc(0.99999,0.00001),z3),ftoc(0.00001,0.0));
   *z2 = csub(cmul(ftoc(1.00001,0.00001),z3),ftoc(0.00001,0.0));
}   
/*=================================================================== */
Zmatrix Unit()
{
  Zmatrix v;
  (v.e11).re = 1.0;
  (v.e11).im = 0.0;
  (v.e21).re = 0.0;
  (v.e21).im = 0.0;
  v.e22 = v.e11;
  v.e12 = v.e21;
  return( v );
}  
/*=================================================================== */
zvector zvm(Zmatrix T, zvector v)
{
  zvector u;

  u.e1 = cadd(cmul(T.e11,v.e1),cmul(T.e12,v.e2));
  u.e2 = cadd(cmul(T.e21,v.e1),cmul(T.e22,v.e2));

  return( u );
}  
