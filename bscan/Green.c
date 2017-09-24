#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complx.h>
#include <jkb.h>
#include "Green.h"
/* File Global variables --------------------------------*/
/* ---Off Diagonal Integration---------------------------*/

static double xi[6]={-0.93246951420315202781,
		     -0.66120938646626451366,
		     -0.23861918608319690863,
		     0.23861918608319690863,
		     0.66120938646626451366,
		     0.93246951420315202781};
static double w[6]={0.17132449237917034504,
		    0.36076157304813860756,
		    0.46791393457269104738,
		    0.46791393457269104738,
		    0.36076157304813860756,
		    0.17132449237917034504};
static int mgpts = 6;

/*
  static double xi[12]={-0.98156063424671925069,
  -0.90411725637047485667,
  -0.76990267419430468703,
  -0.58731795428661744729,
  -0.36783149899818019375,
  -0.12523340851146891547,
  0.12523340851146891547,
  0.36783149899818019375,
  0.58731795428661744729,
  0.76990267419430468703,
  0.90411725637047485667,
  0.98156063424671925069};
  
  static double w[12]={0.047175336386511827194,
  0.10693932599531843096,
  0.16007832854334622633,
  0.20316742672306592174,
  0.23349253653835480876,
  0.24914704581340278500,
  0.24914704581340278500,
  0.23349253653835480876,
  0.20316742672306592174,
  0.16007832854334622633,
  0.10693932599531843096,
  0.047175336386511827194};
  static int mgpts = 12;
*/
/*-----Diagonal Integration------------------------------*/
static double c1 =  0.6366197723675813; /* 2/pi */

static int nns = 6;
static int ns = 6;
static double x[6]={0.12523340851146891547,
		    0.36783149899818019375,
		    0.58731795428661744729,
		    0.76990267419430468703,
		    0.90411725637047485667,
		    0.98156063424671925069};
static double a[6]={0.24914704581340278500,
		    0.23349253653835480876,
		    0.20316742672306592174,
		    0.16007832854334622633,
		    0.10693932599531843096,
		    0.04717533638651182719};
static double xs[6]={0.02163400584411694899,
		     0.12958339115495079613,
		     0.31402044991476550879,
		     0.53865721735180214454,
		     0.75691533737740285216,
		     0.92266885137212023733};
static double as[6]={0.23876366257854756972,
		     0.30828657327394679296,
		     0.24531742656321038598,
		     0.14200875656647668542,
		     0.05545462232488629002,
		     0.01016895869293227589};

/*
  static int ns = 16;
  static int nns  = 16;
  static double x[16] ={0.048307665687738316234,
  0.14447196158279649348,
  0.23928736225213707454,
  0.33186860228212764977,
  0.42135127613063534536,
  0.50689990893222939002,
  0.58771575724076232904,
  0.66304426693021520097,
  0.73218211874028968038,
  0.79448379596794240696,
  0.84966761373256997013,
  0.89632115576605212396,
  0.93490607593773968917,
  0.96476225558750643077,
  0.98561151154526833540,
  0.99726386184948156354};
  static double a[16] ={0.096540088514727800566,
  0.095638720079274859419,
  0.093844399080804565639,
  0.091173878695763884715,
  0.087652093004403811142,
  0.083311924226946755222,
  0.078193895787070306471,
  0.072345494108848506225,
  0.065822222776361846837,
  0.058684093478535547145,
  0.050998059262376176196,
  0.042835898022226680656,
  0.034273862913021433102,
  0.025392065309262059455,
  0.016274394730905670605,
  0.0070186100094700966004};
  static double xs[16]={0.0038978344871159159240,
  0.023028945616873239820,
  0.058280398306240412348,
  0.10867836509105403648,
  0.17260945490984393776,
  0.24793705447057849514,
  0.33209454912991715598,
  0.42218391058194860011,
  0.51508247338146260347,
  0.60755612044772872408,
  0.69637565322821406115,
  0.77843256587326540520,
  0.85085026971539108323,
  0.91108685722227190541,
  0.95702557170354215759,
  0.98704780024798447675};
  static double as[16]={0.060791710043591232851,
  0.10291567751758214438,
  0.12235566204600919355,
  0.12756924693701598871,
  0.12301357460007091542,
  0.11184724485584572262,
  0.096596385152124341252,
  0.079356664351473138782,
  0.061850494581965207095,
  0.045435246507726668628,
  0.031098974751581806409,
  0.019459765927360842078,
  0.010776254963205525645,
  0.0049725428900876417125,
  0.0016782011100511945150,
  0.00028235376466843632177};
*/
/*-----------------------------------------------------*/  
static int na;
static int nc;
static int nb;

static double dd2;
static double dcd2;
static double teta;
static double feta;

static complx g;
static complx *p1;
static complx *q1;
static complx *p2;
static complx *q2;

/* End File Global variables ----------------------------*/

void dodiag(int na, int nc, int n, double etatmp, double del, double delc, 
	    double *g, double *h )
{
  int i, ii;
  double gint;
  
  gint = diag( del, etatmp );
  for( i = 0; i < na; i++ ){
    ii = i*(n+1);
    g[ii] = gint;
    h[ii] = 0.5e0;
  }
  for( i = na + nc; i < (2*na + nc); i++ ){
    ii = i*(n+1);
    g[ii] = gint; 
    h[ii] = 0.5e0;
  }
  gint = diag( delc, etatmp );
  for( i = na; i < (na + nc); i++ ){
    ii = i*(n+1);
    g[ii] = gint;
    h[ii] = 0.5e0;
  }
  for( i = 2*na + nc; i < n; i++ ){
    ii = i*(n+1);
    g[ii] = gint;
    h[ii] = 0.5e0;
  }
  return;
}
/*----------------------------------------------------------------------- */
double diag(double del, double eta)
{
  int i;
  double dd2;
  double gi;
 
  dd2 = 0.5*del;
  gi = 0.e0;
  /* First do the singular integration */
  for( i = 0; i < ns; i++ ){
    gi += as[i]*j0(eta*dd2*xs[i]);
  }
  gi *= -c1;
  /* Now do the nonsingular integration */
  for( i = 0; i < nns; i++ ){
    gi += a[i]*( y0(eta*dd2*x[i]) - c1*j0(eta*dd2*x[i])*log(x[i]) );
  }
  gi *= -0.25e0*del; /* -1/4 comes from def of G */
  return( gi );
}

/*----------------------------------------------------------------------- */
void dorev(int na, int nc, int n1, double *g1, double *h1)
{
  int nl, nr;

  nl = 2*na + nc;
  nr = n1;
  revers( n1, nc, nl, nr, g1, h1 ); /*reverse left side*/
  nl = na + nc;
  nr = 2*na + nc;
  revers( n1, na, nl, nr, g1, h1 ); /*reverse top side*/
  return;
}
/*----------------------------------------------------------------------- */
void revers(int n, int nc, int nl, int nr, double *a, double *b)
{
  int i, ijl, ijr, j, ncd2;
  double tmp;

  ncd2 = nc/2;
  for( j = 0; j < ncd2; j++ ){
    for( i = 0; i < n; i++ ){
      ijl = i + (j + nl)*n;
      ijr = i + (nr - j - 1)*n;
      tmp = a[ijl];
      a[ijl] = a[ijr];
      a[ijr] = tmp;
      tmp = b[ijl];
      b[ijl] = b[ijr];
      b[ijr] = tmp;
    }
  }
  return;
}
/*----------------------------------------------------------------------- */
void green(int na, int nc, int n, int *iz, int *ix, int *ir, 
	   double del, double delc, double eta, double *gi, double *hi)
{
  int lh;
  int i, ij, ijs, j, k;
  double ar, ay0[mgpts], ay1[mgpts], cy1[mgpts], dcd2, dd2, dtmp, rt[mgpts];
  int Have;

  dd2 = 0.5e0*del;
  dcd2 = 0.5e0*delc;
  for( j = 0; j < n; j++ ){
    for( i = 0; i < n; i++ ){
      ij = i + j*n;
      gi[ij] = 0.e0;
      hi[ij] = 0.e0;
      Have = 0;
      if( i != j ){
	for( ijs = 0; ijs < (ij - 1); ijs++ ){
	  if( ir[ij] == ijs ){
	    Have = ijs;
	    break;
	  }
	}
	if( Have ){ /* Changes here */
	  gi[ij] = gi[Have];
	  hi[ij] = hi[Have];
	}
	else{
	  if( j < na ){
	    lh = (ix[j] - ix[i]);
	    for( k = 0; k < mgpts; k++ ){
	      ar = sqrt( SQ(dd2*(iz[j] - iz[i] + xi[k])) + 
			 SQ(dcd2*(ix[j] - ix[i] )) );
	      rt[k] = eta*ar;
	      cy1[k] = eta*dcd2*(ix[j] - ix[i])/ar;
	    }
	    dtmp = dd2;
	  }
	  else if( j < na + nc ){
	    lh = (iz[i] - iz[j]);
	    for( k = 0; k < mgpts; k++ ){
	      ar = sqrt( SQ(dd2*(iz[j] - iz[i])) + 
			 SQ(dcd2*(ix[j] - ix[i] + xi[k])) );
	      rt[k] = eta*ar;
	      cy1[k] = eta*dd2*(iz[i] - iz[j])/ar;
	    }
	    dtmp = dcd2;
	  }
	  else if( j < 2*na + nc ){
	    lh = (ix[i] - ix[j]);
	    for( k = 0; k < mgpts; k++ ){
	      ar = sqrt( SQ(dd2*(iz[j] - iz[i] + xi[k])) + /* c */ 
			 SQ(dcd2*(ix[j] - ix[i])) );
	      rt[k] = eta*ar;
	      cy1[k] = eta*dcd2*(ix[i] - ix[j])/ar;
	    }
	    dtmp = dd2; /* Reverse sign across top */
	  }
	  else if( j < n ){
	    lh = (iz[j] - iz[i]);
	    for( k = 0; k < mgpts; k++ ){
	      ar = sqrt( SQ(dd2*(iz[j] - iz[i])) + 
			 SQ(dcd2*(ix[j] - ix[i] + xi[k])) ); /* c */
	      rt[k] = eta*ar;
	      cy1[k] = eta*dd2*(iz[j] - iz[i])/ar;
	    }
	    dtmp = dcd2; /* Rever sign across left size */
	  }
	  for( k = 0; k < mgpts; k++ ){
	    ay0[k] = y0( rt[k] );
	    gi[ij] += w[k]*ay0[k];
	  }
	  gi[ij] *= -0.25e0*dtmp;
	  if( lh ){
	    for( k = 0; k < mgpts; k++ ){
	      ay1[k] = y1( rt[k] );
	      hi[ij] += w[k]*cy1[k]*ay1[k];
	    }
	    hi[ij] *= -0.25e0*dtmp;
	  }
	}
      }
    }
  }
  return;
} 
/*---------------------------------------------------------------------- */
void boupts(int na, int nb, int nc, int n1, int n2, 
	    int *iz1, int *iz2, int *ix1, int *ix2, int *ir1, int *ir2)
{
  int i, ij, j, js, k, kl, l, ls;

  for( j = 0; j < na; j++ ){
    iz1[j] = 2*j + 1; /* Bottom points */
    ix1[j] = 0;
    iz1[na + nc + j] = 2*na - (2*j + 1);/* Top points */
    ix1[na + nc + j] = 2*nc;
  }
  l = 0;
  for( j = na; j < (na + nc); j++ ){
    l++;
    iz1[j] = 2*na; /* Right Side */
    ix1[j] = 2*l - 1;
    iz1[na + nc + j] = 0; /* Left side */
    ix1[na + nc + j] = 2*nc - (2*l - 1);
  }
  for( j = 0; j < nb; j++ ){
    iz2[j] = 2*na + 2*j + 1;
    ix2[j] = 0;
    iz2[nb + nc + j] = 2*(na + nb) - (2*j + 1);
    ix2[nb + nc + j] = 2*nc;
  }
  l = 0;
  for( j = nb; j < (nb + nc); j++ ){
    ++l;
    iz2[j] = 2*(na + nb);
    ix2[j] = 2*l - 1;
    iz2[nb + nc + j] = 2*na;
    ix2[nb + nc + j] = 2*nc - (2*l - 1);
  }

  /* Now determine which elements repeat */
  /* js = 1 => top, bottom */
  for( j = 0; j < n1; j++ ){
    if( j < na ){
      js = 1;
    }
    else if( j < na + nc ){
      js = 0;
    }
    else if( j < 2*na + nc ){
      js = 1;
    }
    else{
      js = 0;
    }

    for( i = 0; i < n1; i++ ){
      ij = i + j*n1;
      ir1[ij] = ij;
      for( kl = 0; kl < ij; kl++ ){
	l = kl/n1;
	k = kl - l*n1;
	if( abs( ix1[j] - ix1[i] ) == abs( ix1[l] - ix1[k] ) ){
	  if( abs( iz1[j] - iz1[i] ) == abs( iz1[l] - iz1[k] ) ){
	    if( l < na ){
	      ls = 1;
	    }
	    else if( l < na + nc ){
	      ls = 0;
	    }
	    else if( l < 2*na + nc ){
	      ls = 1;
	    }
	    else{
	      ls = 0;
	    }
	    if( js == ls ){
	      ir1[ij]=kl;
	      break;
	    }
	  }
	}
      }
    }
  }
  /* Second group of points */
  for( j = 0; j < n2; j++ ){
    if( j < nb ){
      js = 1;
    }
    else if( j < nb + nc ){
      js = 0;
    }
    else if( j < 2*nb + nc ){
      js = 1;
    }
    else{
      js = 0;
    }
    for( i = 0; i < n2; i++ ){
      ij = i + j*n2;
      ir2[ij] = ij;
      for( kl = 0; kl < ij; kl++ ){
	l = kl/n2;
	k = kl - l*n2;
	if( abs( ix2[j] - ix2[i] ) == abs( ix2[l] - ix2[k] ) ){
	  if( abs( iz2[j] - iz2[i] ) == abs( iz2[l] - iz2[k] ) ){
	    if( l < nb ){
	      ls = 1;
	    }
	    else if( l < nb + nc ){
	      ls = 0;
	    }
	    else if( l < 2*nb + nc ){
	      ls = 1;
	    }
	    else if( l < n2 ){
	      ls = 0;
	    }
	    if( js == ls ){
	      ir2[ij]=kl;
	      break;
	    }
	  }
	}
      }
    }
  }
  return;
}
/*----------------------------------------------------------------------- */
/* Grating Field Functions below.
 * Gn1 is for the tooth region and
 * Gn2 is for the filled region.
 *------------------------------*/
complx Gn1(double x_, double z_)
{
  static complx Z = {0.0,0.0};
  int j, l, l_;
  complx ex, GmH;
  complx sum;
  complx sum1, sum2, sum3, sum4;
  double G, H;
  double xp, xpmx, xpmx2;
  double zp, zpmz, zpmz2;
  double R, Rn;

  sum.re = 0.0;
  sum.im = 0.0;
  /* -----1. Integrate bottom ---------------------*/
  sum1 = Z;
  xp = 0.0;
  xpmx = xp - x_;
  xpmx2= xpmx*xpmx;
  for(l = 0, l_ = 0; l < na ; l++, l_++){
    G = 0.0;
    H = 0.0;
    zp = (2*l+1-na)*dd2 - z_;
    ex = zexp(cfm(g,-zp));
    for( j = 0; j < mgpts; j++ ){
      zpmz = zp + dd2*xi[j];
      zpmz2= zpmz*zpmz;
      R = sqrt(xpmx2 + zpmz2);
      Rn=-xpmx/R;
      G += w[j]*y0(teta*R);
      H +=-w[j]*teta*y1(teta*R)*Rn; /* dy0/dr = - y1 */
    }
    GmH = cmul( csub(cfm(p1[l_],G),cfm(q1[l_],H)) , ex );
    sum1 = cadd( sum1, GmH );
  }
  sum1 = cfm(sum1,-0.25*dd2);
  /* -----2. Right Side ----------------------------*/
  sum2 = Z;
  zp = na*dd2;
  zpmz = zp - z_;
  ex = zexp(cfm(g,-zpmz));
  zpmz2= zpmz*zpmz;
  for(l = 0, l_ = na ; l < nc; l++, l_++){
    G = 0.0;
    H = 0.0;
    xp = (2*l+1)*dcd2 - x_;
    for( j = 0; j < mgpts; j++){
      xpmx = xp+dcd2*xi[j];
      xpmx2 = xpmx*xpmx;
      R = sqrt(xpmx2 + zpmz2);
      Rn= zpmz/R;
      G += w[j]*y0(teta*R);
      H +=-w[j]*teta*Rn*y1(teta*R); /* dy0/dr = - y1 */
    }
    GmH = cmul( csub(cfm(p1[l_],G),cfm(q1[l_],H)), ex );
    sum2 = cadd(sum2,GmH);
  }
  sum2 = cfm(sum2,-0.25*dcd2);
  /* 3. -----Top Side: It is negative! -------------*/
  sum3 = Z;
  xp = 2*dcd2*nc;
  xpmx = xp - x_;
  xpmx2= xpmx*xpmx;
  for( l = 0, l_ = na+nc; l < na; l++, l_++){
    G = 0.0;
    H = 0.0;
    zp = (2*l+1-na)*dd2 - z_;
    ex = zexp(cfm(g,-zp));
    for( j = 0; j < mgpts; j++ ){
      zpmz = zp+dd2*xi[j];
      zpmz2= zpmz*zpmz;
      R = sqrt(xpmx2 + zpmz2);
      Rn= xpmx/R;
      G += w[j]*y0(teta*R);
      H +=-w[j]*teta*y1(teta*R)*Rn; /* dy0/dr = - y1 */
    }
    GmH = cmul(csub(cfm(p1[l_],G),cfm(q1[l_],H)),ex);
    sum3 = cadd(sum3,GmH);
  }
  sum3= cfm(sum3,-0.25*dd2);
  /* -----4. Left side: It is negative!!! ------------*/
  sum4 = Z;
  zp =-na*dd2;
  zpmz = zp - z_;
  ex = zexp(cfm(g,-zpmz));
  zpmz2= zpmz*zpmz;
  for(l = 0, l_ = 2*na+nc ; l < nc; l++, l_++){
    G = 0.0;
    H = 0.0;
    xp = (2*l+1)*dcd2 - x_;
    for( j = 0; j < mgpts; j++){
      xpmx = xp+dcd2*xi[j];
      xpmx2 = xpmx*xpmx;
      R = sqrt(xpmx2 + zpmz2);
      Rn=-zpmz/R;
      G += w[j]*y0(teta*R);
      H +=-w[j]*teta*Rn*y1(teta*R); /* dy0/dr = - y1 */
    }
    GmH = cmul( csub(cfm(p1[l_],G),cfm(q1[l_],H)), ex );
    sum4 = cadd(sum4,GmH);
  }
  sum4 = cfm(sum4,-0.25*dcd2);
  /* ----- Completed 4 regions -----------------------*/
  sum.re = sum1.re+sum2.re+sum3.re+sum4.re;
  sum.im = sum1.im+sum2.im+sum3.im+sum4.im;
  /*  
  printf("\n(%e,%e)  (%e,%e)  (%e,%e)  (%e,%e)\n",
	 sum1.re,sum1.im,sum2.re,sum2.im,sum3.re,sum3.im,sum4.re,sum4.im);
  printf("%f %f--> sum= (%e,%e)\n",x_,z_, sum.re, sum.im);
  */  
  return( sum );
}
/*----------------------------------------------------------------------- */
complx Gn2(double x_, double z_)
{
  int j, l, l_;
  complx sum;
  complx ex, GmH;
  complx sum1, sum2, sum3, sum4;
  double G, H;
  double xp, xpmx, xpmx2;
  double zp, zpmz, zpmz2;
  double R, Rn;

  sum.re = 0.0;
  sum.im = 0.0;
  /* -----1. Integrate bottom ---------------------*/
  sum1.re = 0.0; sum1.im = 0.0;
  xp = 0.0;
  xpmx = xp - x_;
  xpmx2= xpmx*xpmx;
  for(l = 0, l_ = 0; l < nb ; l++, l_++){
    G = 0.0;
    H = 0.0;
    zp = (2*l+1+na)*dd2 - z_;
    ex = zexp(cfm(g,-zp));
    for( j = 0; j < mgpts; j++ ){
      zpmz = zp+dd2*xi[j];
      zpmz2= zpmz*zpmz;
      R = sqrt(xpmx2 + zpmz2);
      Rn=-xpmx/R;
      G += w[j]*y0(feta*R);
      H +=-w[j]*feta*y1(feta*R)*Rn; /* dy0/dr = - y1 */
    }
    GmH = cmul( csub(cfm(p2[l_],G),cfm(q2[l_],H)) , ex );
    sum1 = cadd( sum1, GmH );
  }
  sum1.re *= -0.25*dd2;
  sum1.im *= -0.25*dd2;
  /* -----2. Right Side ----------------------------*/
  sum2.re = 0.0; sum2.im = 0.0;
  zp = (na+2*nb)*dd2;
  zpmz = zp - z_;
  ex = zexp(cfm(g,-zpmz));
  zpmz2= zpmz*zpmz;
  for(l = 0, l_ = nb ; l < nc; l++, l_++){
    G = 0.0;
    H = 0.0;
    xp = (2*l+1)*dcd2 - x_;
    for( j = 0; j < mgpts; j++){
      xpmx = xp+dcd2*xi[j];
      xpmx2 = xpmx*xpmx;
      R = sqrt(xpmx2 + zpmz2);
      Rn= zpmz/R;
      G += w[j]*y0(feta*R);
      H +=-w[j]*feta*y1(feta*R)*Rn; /* dy0/dr = - y1 */
    }
    //pnew = p2[l_];
    GmH = cmul( csub(cfm(p2[l_],G),cfm(q2[l_],H)), ex );
    sum2 = cadd(sum2,GmH);
  }
  sum2.re *= -0.25*dcd2;
  sum2.im *= -0.25*dcd2;
  /* 3. -----Top Side: It is negative! -------------*/
  sum3.re = 0.0; sum3.im = 0.0;
  xp = 2*dcd2*nc;
  xpmx = xp - x_;
  xpmx2= xpmx*xpmx;
  for( l = 0, l_ = nb+nc; l < nb; l++, l_++){
    G = 0.0;
    H = 0.0;
    zp = (2*l+1+na)*dd2 - z_;
    ex = zexp(cfm(g,-zp));
    for( j = 0; j < mgpts; j++ ){
      zpmz = zp+dd2*xi[j];
      zpmz2= zpmz*zpmz;
      R = sqrt(xpmx2 + zpmz2);
      Rn= xpmx/R;
      G += w[j]*y0(feta*R);
      H +=-w[j]*feta*y1(feta*R)*Rn; /* dy0/dr = - y1 */
    }
    GmH = cmul(csub(cfm(p2[l_],G),cfm(q2[l_],H)),ex);
    sum3 = cadd(sum3,GmH);
  }
  sum3.re *= -0.25*dd2;
  sum3.im *= -0.25*dd2;
  /* -----4. Left side: It is negative!!! ------------*/
  sum4.re = 0.0; sum4.im = 0.0;
  zp = na*dd2;
  zpmz = zp - z_;
  ex = zexp(cfm(g,-zpmz));
  zpmz2= zpmz*zpmz;
  for(l = 0, l_ = 2*nb+nc ; l < nc; l++, l_++){
    G = 0.0;
    H = 0.0;
    xp = (2*l+1)*dcd2 - x_;
    for( j = 0; j < mgpts; j++){
      xpmx = xp+dcd2*xi[j];
      xpmx2 = xpmx*xpmx;
      R = sqrt(xpmx2 + zpmz2);
      Rn=-zpmz/R;
      G += w[j]*y0(feta*R);
      H +=-w[j]*feta*y1(feta*R)*Rn; /* dy0/dr = - y1 */
    }
    //pnew = p2[l_];
    GmH = cmul(csub(cfm(p2[l_],G),cfm(q2[l_],H)),ex);
    sum4 = cadd(sum4,GmH);
  }
  sum4.re *= -0.25*dcd2;
  sum4.im *= -0.25*dcd2;
  /* ----- Completed 4 regions -----------------------*/
  sum.re = sum1.re+sum2.re+sum3.re+sum4.re;
  sum.im = sum1.im+sum2.im+sum3.im+sum4.im;

  //printf("(%e,%e)  (%e,%e)  (%e,%e)  (%e,%e)\n",
  // sum1.re,sum1.im,sum2.re,sum2.im,sum3.re,sum3.im,sum4.re,sum4.im);
  return( sum );
}
/*----------------------------------------------------------------------- */
void GratingFields(int Na, int Nc, int Nb,
		   double Del, double Delc,
		   complx *P1, complx *Q1, complx *P2, complx *Q2,
		   int Nx, int Nz, double *xp, double *zp, 
		   complx gamma, complx *Fld,
		   double toothEta, double fillEta)
     /* Data is Stored X first */
{
  double x, z;
  int i, j;
  int ij;

  na = Na;
  nc = Nc;
  nb = Nb;
  dd2 = 0.5*Del;
  dcd2 = 0.5*Delc;
  teta = toothEta;
  feta = fillEta;
  p1 = P1;
  q1 = Q1;
  p2 = P2;
  q2 = Q2;
  g = gamma;

  for( j = 0; j < Nz; j++){
    for(i = 0; i < Nx; i++){
      ij = i + j*Nx;
      x = xp[i];
      z = zp[j];
      if( j < Na )
	Fld[ij] = Gn1(x,z);
      else
	Fld[ij] = Gn2(x,z);
    }
  }
  return;
}

/*---------------------------------------------------------------*/
void GratingSketch(int nitter, int nflt, int na, int nb, int nc, 
		   double delc, double del, double fko,
		   int GLbot, int GLtop, double *xpts,
		   double *xgrat, double *zgrat, complx *Fgrat,
		   complx g, double ng, double nfill,
		   complx *Vtop, complx *Vbot,
		   complx *p1, complx *q1, complx *p2, complx *q2 )
     /*----------------------------------------------------------
      * Input parameters include:
      * nitter - The itteration number
      * nflt - number of space harmonics
      * na - number of panels at bottom and top of tooth
      * nb - number of panels at bottom and top of fill region
      * nc - number of panels along the sides
      * delc - panel heights
      * del - panel widths
      * fko - free-space wavenumber
      * xgrat - Array of points along x (center of panels)
      * zgrat - Array of points along z (center of panels)
      * Fgrat - Array of complex fields at panel intersections
      * g - complex propagation constant
      * ng - refractive index of tooth
      * nfill - refractive index of fill region
      * Vbot - panel fields at bottom of grating
      * Vtop - panel fields at top of grating
      * p1 - Array of flux values around the tooth retion
      * q1 - Array of field values around the tooth region
      * p2 - Array of flux values around the fill region
      * q2 - Array of field values around the fill region
      *----------------------------------------------------------*/
{
  char FIELD[64], INT[4];
  int i, il, l, n;

  double Maggrat, Phagrat;

  FILE *fMag;
  FILE *fPhase;


  strcpy(FIELD, "Grating/Mag");
  sprintf(INT,"%3d",nitter);
  for(n = 0 ; INT[n] != '\0' ; n++)
    if( INT[n] == ' ')
      INT[n] = '0';
  strcat(FIELD,INT);
  strcat(FIELD,".prn");
  fMag = fopen(FIELD,"w");
  if(fMag == NULL){
    fprintf(stderr, "Could not open '%s' in GratingSketch.\n",FIELD);
    exit( EXIT_FAILURE );
  }

  strcpy(FIELD, "Grating/Pha");
  strcat(FIELD,INT);
  strcat(FIELD,".prn");
  fPhase = fopen(FIELD,"w");
  if(fPhase == NULL){
    fprintf(stderr, "Could not open '%s' in GratingSketch.\n",FIELD);
    exit( EXIT_FAILURE );
  }

  for( l = 0; l < nflt; l++)
    zgrat[l] = (2*l+1-na)*del*0.5; /* Center of tooth!!! */
  for( i = 0; i < nc ; i++)
    xgrat[i] = (2*i+1)*delc*0.5;
  GratingFields(na, nc, nb, del, delc,p1, q1, p2, q2,
		nc, nflt, xgrat, zgrat, 
		g, Fgrat, ng, nfill);

  for( l = 0; l < nflt ; l++){
    Maggrat = zabs2( Vbot[l] );
    Phagrat = atan2( Vbot[l].im,Vbot[l].re );
    fprintf( fMag,"%f\t%f\t%f\n", xpts[GLbot]/fko, zgrat[l]/fko, Maggrat );
    fprintf( fPhase,"%f\t%f\t%f\n", xpts[GLbot]/fko, zgrat[l]/fko, Phagrat );
    for( i = 0; i < nc; i++){
      il = i + nc*l;
      Maggrat = zabs2( Fgrat[il] );
      Phagrat = atan2( Fgrat[il].im,Fgrat[il].re );
      //printf("%f\t%f\t%f\n", xgrat[i], zgrat[l], Maggrat[il]);
      fprintf(fMag,"%f\t%f\t%f\n", xgrat[i]/fko, zgrat[l]/fko, Maggrat );
      fprintf(fPhase,"%f\t%f\t%f\n", xgrat[i]/fko, zgrat[l]/fko, Phagrat );
    }
    Maggrat = zabs2( Vtop[l] );
    Phagrat = atan2( Vtop[l].im,Vtop[l].re );
    fprintf( fMag, "%f\t%f\t%f\n", xpts[GLtop-1]/fko, zgrat[l]/fko, Maggrat );
    fprintf( fPhase, "%f\t%f\t%f\n", xpts[GLtop-1]/fko, zgrat[l]/fko, Phagrat );
    //printf("\n");
    fprintf(fMag,"\n");
    fprintf(fPhase,"\n");
  }
  fclose(fMag);
  fclose(fPhase);
}







