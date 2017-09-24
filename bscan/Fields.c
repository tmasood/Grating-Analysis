/*----------------------------------------------------------------------
 *     Copyright Jerome K. Butler
 *     1990-2001                                        November 2001 
 *----------------------------------------------------------------------*/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complx.h>
#include <jkb.h>
#include "Dft.h"
#include "Green.h"

/*-------------Field.c Globals-----------------------------------*/
static int NL;        /* number of layers */
static int NF;        /* number of space harmonics */
static int Nfl;       /* Space harmonics (-nfl,...,nfl-1)*/
static int Ntop;      /* first or top grating layer */
static int Nbot;      /* last or bottom grating layer */
static int nspl;      /* number of spline points in grating layer (nc+2) */

static double tauz;   /* 2 \tau_z = Normalized toothwidth */

static double *xpts;  /* Layer boundaries */
static double *xrpts; /* Reverse layer boundary numbers */
static double *zpts;  /* zpts[l] = z_l, l = 0 \ldots nflt */
static double *d;     /* Nommalized layer thicknesses */
static double *Gamma; /* Layer confinement factor */
static complx *W;     /* W[n] = exp(-j 2*Pi*n/N) */
static complx *PF;    /* Phase Factor W_{1/2}^{-n(N_a-1)} */
static complx *h;     /* Lateral space harmonics NF x NL dim */
static complx *c;     /* cos(hn*d) NF x NL dim */
static complx *s;     /* sin(hn*d) NF x NL dim */
static complx *t;     /* tan(hn*d) NF x NL dim */
static complx *psi;   /* psin (wave functions) NF x NL dim */
static complx *psip;  /* psinp (derivative( NF x NL dim */
/* Arrays related to the center of tooth fields */
static double *xspl;
static double *psre;
static double *psim;
static double *bre;
static double *cre;
static double *dre;
static double *bim;
static double *cim;
static double *dim;

/*------------------------------------------------------------*/
complx Psi(double );
double Psi2(double );
void Psil(double , complx *, complx *);

/*------------------------------------------------------------*/
double Norm(int nlayers, int nfl, int nflt, int ntop, int nbot,
	    char *comments,
	    int *TopLeak, int *BotLeak,
	    double *dum, double *Xpts, double *Xrpts, 
	    double *zpts_, double tauz_, double *Gam,
	    complx *ht, complx *ct, complx *st, complx *tt,
	    complx *ps, complx *psp, complx *W_, complx *PF_,
	    int nspl_,
	    double *xspl_, double *psre_, double *psim_,
	    double *bre_, double *cre_, double *dre_,
	    double *bim_, double *cim_, double *dim_)
{
  int i, l, m, ml, m_, n, nl, n_;
  double srtotal, total, xhigh, xlow;
  complx F;

  static double eps = 1.e-6;
  static double delta = 1.e-12;

  /* Transfer the variables to make global in fields (see above). */
  NL = nlayers;
  NF = nflt;
  Nfl = nfl;
  Ntop = ntop;
  Nbot = nbot;
  d = dum;
  xpts = Xpts;
  xrpts = Xrpts;
  zpts = zpts_;
  tauz = tauz_;
  Gamma = Gam;
  h = ht;
  c = ct;
  s = st;
  t = tt;
  psi = ps;
  psip = psp;
  W = W_;
  PF = PF_;
  nspl = nspl_;
  xspl = xspl_;
  psre = psre_;
  psim = psim_;
  bre = bre_;
  cre = cre_;
  dre = dre_;
  bim = bim_;
  cim = cim_;
  dim = dim_;
  /* End Global transfer */

  /* confinement factors */
  /*------------------------------------------*/
  /* Bottom layer NL-1 */
  l = NL-1;
  Gamma[l] = 0.0;
  for( n = -nfl; n < nfl; n++){
    n_ = n + nfl;
    nl = n_;
    if( !BotLeak[n_] )
      Gamma[l] += 0.5*(SQ(psi[nl].re)+SQ(psi[nl].im))/h[nl].re;
  }
  for( n = -nfl; n < nfl; n++){
    n_ = n + nfl;
    nl = n_;
    for(m = n+1; m < nfl ; m++){
      m_ = m + nfl;
      ml = m_;
      if( !BotLeak[n_] && !BotLeak[m_]){
	F = tdiv(psi[nl],dconjg(psi[ml]),cadd(h[nl],dconjg(h[ml])));
	Gamma[l] += 2.0*F.re;
      }
    }
  }
   /* Layers below the grating */
  for( l = NL-2 ; l > Nbot ; l-- ) {
    xlow = xpts[l]+delta;
    xhigh = xpts[l-1]-delta;
    Gamma[l] = simp_( Psi2, xlow, xhigh, eps );
  }
  /* Grating layers */
   for( l = Nbot ; l > Ntop-1 ; l-- ) {
    xlow = xpts[l]+delta;
    xhigh = xpts[l-1]-delta;
    Gamma[l] = simp_( Psi2, xlow, xhigh, eps );
  }
 /* Layers above the grating */
  for( l = Ntop-1; l > 0 ; l-- ) {
    xlow = xpts[l]+delta;
    xhigh = xpts[l-1]-delta;
    Gamma[l] = simp_( Psi2, xlow, xhigh, eps );
  }
  /* Top layer 0 */
  l = 0;
  Gamma[l] = 0.0;
  for( n = -nfl; n < nfl; n++) {
    n_ = n + nfl;
    nl = n_;
    if( !TopLeak[n_] )
      Gamma[l] += 0.5*(SQ(psi[nl].re)+SQ(psi[nl].im))/h[nl].re;
  }
  for( n = -nfl; n < nfl; n++) {
    n_ = n + nfl;
    nl = n_;
    for(m = n+1; m < nfl ; m++) {
      m_ = m + nfl;
      ml = m_;
      if( !TopLeak[n_] && !TopLeak[m_]){
	F = tdiv(psi[nl],dconjg(psi[ml]),cadd(h[nl],dconjg(h[ml])));
	Gamma[l] += 2.0*F.re;
      }
    }
  }
  /*------------------------------------------------------------*/    
  total = 0.0;
  for( l = 0; l < nlayers; l++ )
    total += Gamma[l];
  
  srtotal = 1.0/sqrt(total);
  for( l = 0; l < nlayers; l++ )
    Gamma[l] /= total;

  /* Normalize Space harmonic amplitudes */
  for(l = 0; l < nlayers; l++)  
    for(n = 0; n < NF ; n++){
      nl = n + NF*l;
      ps[nl].re *= srtotal;
      ps[nl].im *= srtotal;
      psp[nl].re *= srtotal;
      psp[nl].im *= srtotal;
    }

  /* Normalize Grating Splined fields */
  for(i = 0; i < nspl; i++) {
    psre[i] *= srtotal;
    psim[i] *= srtotal;
    bre[i] *= srtotal;
    cre[i] *= srtotal;
    dre[i] *= srtotal;
    bim[i] *= srtotal;
    cim[i] *= srtotal;
    dim[i] *= srtotal;
  }
   
  printf( "\n Layer Confinement Factors \n");
  for( l = 0; l < nlayers ; ++l)
    printf( "Gamma(%2d) = %15.4e %s\n",l,Gamma[l],&comments[128*l]);

  /*     check the integral of the fourier transform
   *
   *      xlow =-1.d0
   *      xhigh= 1.d0
   *      call simp_(far,xlow,xhigh,1.d-2,total)
   *      write(6,*)'far-field sum/2pi',total/tpi */
  return ( srtotal );
} 

/*----------------------------------------------------------------
 *
 *  This function evaluates the field intensity at z = 0;
 *
 *    Psi2 = psi*psi
 *  where
 *    psi  = Sum psi_n(u).
 *
 *  The boundary points along u (the x variable) are located at
 *  x(i) where
 *      x(i) .lt. u .lt. x(i+1), (x = xrpts)
 *
 *  input..
 *
 *    n = layers-1, the number of data points
 *    u = the abscissa at which the intensity is to be evaluated
 *
 *  if  u  is not in the same interval as the previous call, then a
 *  binary search is performed to determine the proper interval.
 *------------------------------------------------------------*/
complx Psi(double u)
{
  int j, k, l, n, n_, nl;
  static int i = 0;
  double dx;
  complx field, ft, hnx, ct, tt, zcn, zsn;

  field.re = 0.0;
  field.im = 0.0;

  if(u < xrpts[0] ) { /* Substrate ----------------------*/
    i = -1;
    l = NL - i - 2;
    dx = u - xrpts[0];
    for(n = -Nfl; n < Nfl; n++) {
      n_ = Nfl + n;
      nl = n_ + l*NF;
      hnx = fcm(dx,h[nl]);
      field = cadd(field,cmul(psi[nl],zexp(hnx)));
    }
    return( field );
  } /* End Substrate */

  if( u < xrpts[NL-Nbot-2] ) { /* Layers below grating */
    if( i >= NL - 1 )
      i = 0;
    if( (u < xrpts[i]) || (u >= xrpts[i+1]) ){      
      i = 0;
      j = NL-2;
      while( j != i+1 ){ 
	k = (i + j)/2;
	if( u < xrpts[k] )
	  j = k;
	else{   
	  if( u >= xrpts[k] )
	    i = k;
	}
      }
    }
    l = NL - i - 2;
    dx = u - xrpts[i];
    for( n = -Nfl; n < Nfl; n++) {
      n_ = Nfl + n;
      nl = n_ + l*NF;
      hnx = fcm(dx,h[nl]);
      ct = ccos(hnx);
      tt = cdiv(ztan(hnx),h[nl]);
      zcn = psi[nl];
      zsn = cmul(tt,psip[nl]);
      ft = cmul(ct,cadd(zcn,zsn));
      field.re += ft.re;
      field.im += ft.im;
    }
    return( field );
  } /* End for Layers below grating -----------------------*/

  if( u < xrpts[NL-Ntop-1] ) { /* Grating Layer -----------*/
    l = Nbot;
    i = NL - l - 2;
    /* field.re = seval(nspl, u, xspl, psre, bre, cre, dre);
    field.im = seval(nspl, u, xspl, psim, bim, cim, dim); */
    return ( field );
  } /* End grating Layer ---------------------------------*/

  if( u < xrpts[NL-2] ) { /* Layers above Grating Layer --*/
    if( i >= NL - 1 )
      i = 0;
    if( (u < xrpts[i]) || (u >= xrpts[i+1]) ) {      
      i = 0;
      j = NL-2;
      while( j-i != 1 ){ 
	k = (i + j)/2;
	if( u < xrpts[k] )
	  j = k;
	else{   
	  if( u >= xrpts[k] )
	    i = k;
	}
      }
    }   
    l = NL - i - 2;
    dx = u - xrpts[i];
    for( n = -Nfl; n < Nfl; n++) {
      n_ = Nfl + n;
      nl = n_ + l*NF;
      hnx = fcm(dx,h[nl]);
      ct = ccos(hnx);
      tt = cdiv(ztan(hnx),h[nl]);
      zcn = psi[nl];
      zsn = cmul(tt,psip[nl]);
      ft = cmul(ct,cadd(zcn,zsn));
      field.re += ft.re;
      field.im += ft.im;
    }
    return( field );
  } /* End for Layers above grating */

  else{ /* Superstrate region */
    l = 0;
    i = NL - l - 2;
    dx = u - xrpts[i];
    for(n = -Nfl; n < Nfl; n++) {
      n_ = Nfl + n;
      nl = n_;
      hnx = fcm(-dx,h[nl]);
      field = cadd(field,cmul(psi[nl],zexp(hnx)));
    }
    return(field);
  } /* End Superstrate region*/
  
  return( ftoc(0.0,0.0) );
}
/*---------------------------------------------------------------*/
double Psi2(double x)
{
  double f;
  complx F;

  F = Psi(x);
  f = F.re*F.re+F.im*F.im;

  //  printf(" %e  %e\n", x, f);
  return (f);
}

/*---------------------------------------------------------------*/
void ModeSketch(int nitter, int N, double *X, double ko,
		complx *vn, complx *Vl)  
     /*---------------------------------------
      * out - the directory i.e. Grating 
      * nitter - the itteration number 
      * N - Number of X points
      * X - Array of x points 
      * ko - free-space wavelength
      * vn - Array working vector (space harmonics)
      * VL - Array working vector ( field at z points )
      *---------------------------------------*/
{
  char FIELD[64], INT[4];

  int i, l, n;
  double x;
  complx F;

  FILE *fields;

  strcpy(FIELD, "Loaded/");
  strcat(FIELD,"field");
  sprintf(INT,"%3d",nitter);
  for(n = 0 ; INT[n] != '\0' ; n++)
    if( INT[n] == ' ')
      INT[n] = '0';
  strcat(FIELD,INT);
  strcat(FIELD,".prn");
  fields = fopen(FIELD,"w");
  /* 1D Field at z = 0 */
  if(nitter < 1000){
    for (i = 0; i < N; i++){
      x = X[i];
      F = Psi(x);
      fprintf(fields," %f\t%f\n",x/ko, zabs2( F ) );
      //      printf("%f  %e\n", x/ko, zabs(F));
    }
  }
  fclose(fields);
  /* 2D fields for 1 period \Clambda */
  strcpy(FIELD, "Loaded/");
  strcat(FIELD,"2D");
  sprintf(INT,"%3d",nitter);
  for(n = 0 ; INT[n] != '\0' ; n++)
    if( INT[n] == ' ')
      INT[n] = '0';
  strcat(FIELD,INT);
  strcat(FIELD,".prn");
  fields = fopen(FIELD,"w");

  if(nitter < 1000){
    for (i = 0; i < N; i++){
      x = X[i];
      Psil(x, vn, Vl);
      for( l = 0; l < NF; l++) {
	fprintf(fields," %f\t%f\t%f\n",x/ko, zpts[l]/ko, zabs( Vl[l] ) );
      }
      fprintf(fields, "\n");
    }
  }
  fclose(fields);
}
/*----------------------------------------------------------------
 *
 *  Input:  psin  = psi_n(u)*PF_n
 *  Output: psil = idft(psin)
 *
 *  where the phase factor is
 *    W_{1/2}^{-n(N_a-1)
 *
 *  The boundary points along u (the x variable) are located at
 *  x(i) where
 *      x(i) .lt. u .lt. x(i+1), (x = xrpts)
 *
 *  input..
 *
 *    n = layers-1, the number of data points
 *    u = the abscissa at which the intensity is to be evaluated
 *
 *  if  u  is not in the same interval as the previous call, then a
 *  binary search is performed to determine the proper interval.
 *------------------------------------------------------------*/
void Psil(double u, complx *psin, complx *psil )
{
  int j, k, l, ll, n, n_, nl;
  static int i = 0;
  double dx, ztemp;
  complx field, ft, hnx, ct, tt, zcn, zsn;

  field.re = 0.0;
  field.im = 0.0;

  if(u < xrpts[0] ) { /* Substrate ----------------------*/
    i = -1;
    l = NL - i - 2;
    dx = u - xrpts[0];
    for(n = -Nfl; n < Nfl; n++) {
      n_ = Nfl + n;
      nl = n_ + l*NF;
      hnx = fcm(dx,h[nl]);
      psin[n_] = tmul( PF[n_], psi[nl], zexp(hnx) );
    }
    dft(NF, Nfl, W, psin, psil);
    return;
  } /* End Substrate */

  if( u < xrpts[NL-Nbot-2] ) { /* Layers below grating */
    if( i >= NL - 1 )
      i = 0;
    if( (u < xrpts[i]) || (u >= xrpts[i+1]) ){      
      i = 0;
      j = NL-2;
      while( j != i+1 ){ 
	k = (i + j)/2;
	if( u < xrpts[k] )
	  j = k;
	else{   
	  if( u >= xrpts[k] )
	    i = k;
	}
      }
    }
    l = NL - i - 2;
    dx = u - xrpts[i];
    for( n = -Nfl; n < Nfl; n++) {
      n_ = Nfl + n;
      nl = n_ + l*NF;
      hnx = fcm(dx,h[nl]);
      ct = ccos(hnx);
      tt = cdiv(ztan(hnx),h[nl]);
      zcn = psi[nl];
      zsn = cmul(tt,psip[nl]);
      ft = cmul(ct,cadd(zcn,zsn));
      psin[n_] = cmul(PF[n_],ft);
    }
    dft(NF, Nfl, W, psin, psil);
    return;
  } /* End for Layers below grating -----------------------*/

  if( u < xrpts[NL-Ntop-1] ) { /* Grating Layer -----------*/
    l = Nbot;
    i = NL - l - 2;
    for( ll = 0; ll < NF; ll++) {
      ztemp = zpts[ll];
      if( ztemp < tauz )
	psil[ll] = Gn1(u,ztemp);
      else
	psil[ll] = Gn2(u,ztemp);
    }
    return;
  } /* End grating Layer ---------------------------------*/

  if( u < xrpts[NL-2] ) { /* Layers above Grating Layer --*/
    if( i >= NL - 1 )
      i = 0;
    if( (u < xrpts[i]) || (u >= xrpts[i+1]) ) {      
      i = 0;
      j = NL-2;
      while( j-i != 1 ){ 
	k = (i + j)/2;
	if( u < xrpts[k] )
	  j = k;
	else{   
	  if( u >= xrpts[k] )
	    i = k;
	}
      }
    }   
    l = NL - i - 2;
    dx = u - xrpts[i];
    for( n = -Nfl; n < Nfl; n++) {
      n_ = Nfl + n;
      nl = n_ + l*NF;
      hnx = fcm(dx,h[nl]);
      ct = ccos(hnx);
      tt = cdiv(ztan(hnx),h[nl]);
      zcn = psi[nl];
      zsn = cmul(tt,psip[nl]);
      ft = cmul(ct,cadd(zcn,zsn));
      psin[n_] = cmul(PF[n_],ft);
    }
    dft(NF, Nfl, W, psin, psil);
    return;
  } /* End for Layers above grating */

  else{ /* Superstrate region */
    l = 0;
    i = NL - l - 2;
    dx = u - xrpts[i];
    for(n = -Nfl; n < Nfl; n++) {
      n_ = Nfl + n;
      nl = n_;
      hnx = fcm(-dx,h[nl]);
      psin[n_] = tmul(PF[n_], psi[nl], zexp(hnx) );
    }
    dft(NF, Nfl, W, psin, psil);
    return;
  } /* End Superstrate region*/
}
