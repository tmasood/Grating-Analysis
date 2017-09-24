/*     Copyright Jerome K. Butler
 *     1990-2001                                           September 2001 
 *-----------------------------------------------------------------------
 *     n(i) refractive index of            ith layer
 *     k(i) extinction coefficient of      ith layer
 *     alpha(i) wave absorption /cm of     ith layer
 *        lossy material entered as a negative real.
 *-----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <jkb.h>
#include <complx.h>
/*----------------Start-Function Prototypes------------------------------*/
void coef(int , 
	  complx *, complx *,
	  complx *, complx *, complx *);

void norm(double (* )(double ), int , 
	  complx , complx *, complx *, complx *);

void points(int , int , double *, double *, double *);
void testdet(double , double , double , double );
void IndexProfile(int , double , double *, complx *);

complx Phih(double );
double  Phih2(double );

complx	fh(complx ), fq(complx ),
	zout(complx , complx , complx , int ,
	     char *,
	     complx *, complx *,
	     complx *, complx *, complx *,
	     double (* )(double ),
	     complx (* )(double ) );

Zmatrix Zmul(Zmatrix , Zmatrix ),
	Ztrans(complx , double ),
	Ztransform(complx , double );

/*------------------End-Function Prototypes------------------------------*/
/*------------------File Global Variables--------------------------------*/
char *comments; /* Layer comments for information */

int NL;
int GLBottom; /* Bottom of the grating */
int GLTop;    /* Top of the grating */

double ko; /* Free-space wavenumber */
double *d; /* Normalized (to ko) layer thicknesses */
double *xpts;
double *xrpts;
double *Xf;
double *gammaz;

complx *kappaz;
complx *hz;
complx *czh;
complx *szh;
complx *phih;
complx *phiph;

/* static complx ZZERO = {0.0, 0.0};
static complx ZONE = {1.0, 0.0};
static complx ZJ = {0.0, 1.0}; */

/*-----------------------------------------------------------------------*/
complx znl(char *comments_,
	   int itmax, 
	   int NL_, 
	   int *GL, 
	   int *DC,
	   double ko_,
	   double *d_,
	   double *xrpts_,
	   double *xpts_,
	   double *Xf_,
	   double *gammaz_,
	   complx *kappaz_,
	   complx *hz_,
	   complx *czh_,
	   complx *szh_,
	   complx *phi_,
	   complx *phip_,
	   complx zg)
/***********************************************************************
 * itmax is the allowed number of iterations in Mueller.
 * Scan and CScan are used to open files for output such as the
 * field distribution.
 **********************************************************************/
{
   char FIELD[32];
   int itter, l, root;
   complx z[3], zw, zroot;

   /* set up the array pointers */
   comments = comments_;
   NL = NL_;
   d = d_;
   xpts = xpts_;
   xrpts = xrpts_;
   Xf = Xf_;
   gammaz = gammaz_;
   ko = ko_;
   kappaz = kappaz_;
   hz = hz_;
   czh = czh_;
   szh = szh_;
   phih = phi_;
   phiph= phip_;

  for(l = 0; l < NL; l++){
    GLTop = l;
    if(GL[l])
      break;
  }
  for(l = NL-2; l > 0; l--){
    GLBottom = l;
    if(GL[l])
      break;
  }

   points(GLBottom,NL,d,xpts,xrpts);
   IndexProfile(NL, ko, xpts, kappaz);

   fprintf(stdout, "*************************************************************************\n");
   fprintf(stdout, "*            Results for Mode With Grating Unloaded\n");
   fprintf(stdout, "***********************************w(z)=0********************************\n");
   fprintf(stdout, "      Real(z)               Imag(z)           Real(w)        Imag(w)\n");
   *(z+2) = zg;
   cguesses(z);
   zroot = croot( fh, z, &zw, &itter, itmax, &root ); /* Replace zg with the root. */

   if( root ){
     strcpy(FIELD,"./");
     strcat(FIELD,"Unloaded");
     strcat(FIELD,"/");
     zout(zg, *(z+2), zw, itter,FIELD, czh, szh, hz, phih, phiph, Phih2, Phih);
     }
   else{
      fprintf( stdout, "Cannot find root with present guess.\n" );
      exit( EXIT_FAILURE );
      }
      
   return( zroot );
}
/*------------------------------------------------------------------------*/
complx zout( complx zg, complx zr, complx fvalue, int itt,
              char *Scanmode,
              complx *CZ, complx *SZ,
              complx *H, complx *PS, complx *PSP,
              double (*fncn)(double ),
              complx (*zfncn)(double) )
{
   int i, l;
   double deltax,x, xstart, xend;
   double TWO, ystart, ystart0, yend, yend0;
   complx F;
   FILE *fields;

   static int NFld = 100;

   fprintf( stdout, "\nConverged in %i Iterations.\n",itt);
   fprintf( stdout, " zguess   = ( %15.10e , %15.10f ) \n",zg.re,zg.im);
   fprintf( stdout, " zroot    = ( %20.15e , %20.15f ) \n",zr.re,zr.im);
   fprintf( stdout, " w(zroot) = ( %15.10e , %15.10e ) \n",fvalue.re,fvalue.im);
   coef(NL, CZ, SZ, H, PS, PSP);
   norm(fncn, NL, zr, H, PS, PSP);

   fprintf(stdout, "\n");

   for( i = 0 ; i < NL-1 ; ++i)
      fprintf(stdout, "phi (%2d) = %12.5e,%12.5e   phip(%2d) = %12.5e,%12.5e \n",
                 i,PS[i].re,PS[i].im,i,PSP[i].re,PSP[i].im );
/*
  printf("h values\n");
  cout_(NL,1,H);               

  printf("czh values\n");
  cout_(NL,1,CZ);               

  printf("szh values\n");
  cout_(NL,1,SZ);               
*/
  strcat(Scanmode,"field.prn");
  fields = fopen(Scanmode,"w");
  if(fields == NULL){
    fprintf(stderr, "Could not open '%s'\n",Scanmode);
    exit(0);
  }
  /* Fix the start and end amplitutes the same. This requires
   * the manipulation of the amplitudes at the beginning and
   * end of the start and end layers.                       */
  TWO = 1.0;
  ystart0 = zabs(PS[NL-2]);
  yend0 = zabs(PS[0]);
  if( ystart0 > yend0 ) {
    xend = xpts[0] + TWO/(H[0].re);
    yend = yend0*exp(-TWO);
    xstart = xpts[NL-2] + log(yend/ystart0)/(H[NL-1].re);
  }
  else {
    xstart = xpts[NL-2]-TWO/H[NL-1].re;  /* Minus sign => below grating */
    ystart = ystart0*exp(-TWO);
    xend = xpts[0] + log(ystart/yend0)/H[0].re;
  }
  deltax = (xend - xstart)/NFld;
  x = xstart-deltax;
  for (l = 0; l < NFld; l++) {
    x += deltax;
    Xf[l] = x;
    F = zfncn(x);
    //fprintf(fields," %f   %f %f\n",x/ko,F.re,F.im);
    fprintf(fields," %f\t%f\n", x/ko, zabs2(F) );
  }
  fclose(fields);
  return( zr );
}
/*-----------------------------------------------------------------------*/
complx fh(complx g)
{
  int l;
  complx arg, c, s, g2, e1, e2, v, ht2, hb2, t, v1, v2;

  g2.re = g.re*g.re - g.im*g.im;
  g2.im = 2.0*g.re*g.im;

  ht2 = cneg(cadd(g2,kappaz[0]));
  hb2 = cneg(cadd(g2,kappaz[NL-1]));

  hz[0] = zsqrt( ht2 );
  hz[NL-1] = zsqrt( hb2 );

  for( l = 1; l < (NL-1); l++ )
    hz[l] = zsqrt( cadd(kappaz[l],g2) );

  e1 = ZONE;
  e2 = *(hz+NL-1);
  
  for(l = NL-2; l > 0; l--){
    arg = fcm(*(d+l),*(hz+l));
    if( fabs(arg.im) < 1.0 ){
      c = ccos(arg);
      s = csin(arg);
      v1 = cadd(cmul(c,e1),tdiv(e2,s,*(hz+l)));
      v2 = csub(cmul(c,e2),tmul(e1,s,*(hz+l)));
      }
    else{
      t = ztan(arg);
      v1 = cadd(e1,tdiv(e2,t,*(hz+l)));
      v2 = csub(e2,tmul(e1,t,*(hz+l)));
      }
    e1 = v1;
    e2 = v2;
    }

  v = cadd(e2,cmul(*hz,e1));
  
  printf("%20.15f  %18.10E  %13.5E  %13.5E \n",g.re,g.im,v.re,v.im);
  return( v );
} 

/*-----------------------------------------------------------------------*/
Zmatrix Ztransform(complx h, double d)
{
   Zmatrix v;
   complx hd, t;

   hd.re = h.re*d;
   hd.im = h.im*d;

   if( fabs(hd.im) < 10.0){
      t = ztan( hd );
      v.e11 = ZONE;
      v.e21 = cneg(cmul(t,h));
      v.e12 = cdiv(t,h);
      v.e22 = ZONE;
      }
   else {
      t = btan(hd);
      v.e11 = ZZERO;
      v.e21 = cneg(cmul(t,h));
      v.e12 = cdiv(t,h);
      v.e22 = ZZERO;
      }
      
   return(v);
} 
/*------------------------------------------------------------------------*/
/************************************************************************
 *     The function calcuates the field coefficients for
 *     the field in the lth layer is given by:
 *
 *     phi(x) = phi(l)*cos[h(x-x(l))] + phip(l)*sin[(h(x-xlay)]/h
 ************************************************************************ */
void coef(int L,
          complx *cz, complx *sz,
          complx *H, complx *ri, complx *rip)
{
  int l;
  double fmag, fmax;
  complx Fmax;
  Zmatrix Oz;

  /*     calculate the vector phi(i) = col{phi(i),phip(i)}
   * */
   
  cz[0] = ZZERO;
  sz[0] = ZZERO;
  for(l = 1; l < L -1 ; l++){
    cz[l] = ccos(fcm(d[l],H[l]));
    sz[l] = csin(fcm(d[l],H[l]));
    }
  cz[L-1] = ZZERO; 
  sz[L-1] = ZZERO; 

  ri[L-1] = ZONE;
  rip[L-1] = H[L-1];
  ri[L-2] = ri[L-1];
  rip[L-2] = rip[L-1];
  for( l = L-2 ; l > 0; --l ){
    Oz = Ztransform(H[l],d[l]);
    ri[l-1] = cadd(tmul(cz[l],Oz.e11,ri[l]),tmul(cz[l],Oz.e12,rip[l]));
    rip[l-1]= cadd(tmul(cz[l],Oz.e21,ri[l]),tmul(cz[l],Oz.e22,rip[l]));
    }

  /* Find the maximum field at the boundaries and normalize relative
     to the maximum. This will be a first step toward the normalized
     values which will be done by integration.*/  
  fmax = 0.;
  Fmax.re = 0., Fmax.im = 0.;
  for( l = 0; l < L ; l++){
    fmag = zabs(ri[l]);
    if( fmag > fmax ){
      fmax = fmag;
      Fmax = ri[l];
      }
    } 
    
  for( l = 0; l < L; l++){
    ri[l] = cdiv(ri[l],Fmax);
    rip[l]= cdiv(rip[l],Fmax);
    }   

  fprintf(stdout, "\n Seminormalized fields at boundaries.\n");
  for( l = 0; l < L-1 ; ++l)
    fprintf(stdout, " phi(%2d) = %15.8e %15.8e \n",
                         l, (ri[l]).re, (ri[l]).im);

  return;
}/*END coef*/ 
/*--------------------------------------------------------------------------*/
void norm(double (*fncn)(double ), int L, 
          complx gamma, complx *H, complx *ri, complx *rip)
{
   int l;
   double delta, eps, srtotal, total, xhigh, xlow;

   eps = 1.e-8;
   delta = 1.e-12;

   l = 0;
   if( H[l].re < 0 ){
      fprintf(stdout, "\n************************************************\n");
      fprintf(stdout,   "* Mode leaks in the top layer. NOT NORMALIZED! *\n");
      fprintf(stdout,   "************************************************\n\n");
      fprintf(stdout, "Leaks at %f Degrees\n", 180.0*atan(gamma.im/H[l].im)/M_PI);
      gammaz[l] = 0.;
      }
   else
      gammaz[l] = 0.5*(ri[l].re*ri[l].re+ri[l].im*ri[l].im)/H[l].re;

   l = L - 1;
   if( H[l].re < 0 ){
      fprintf(stdout, "\n***************************************************\n");
      fprintf(stdout,   "* Mode leaks in the bottom layer. NOT NORMALIZED! *\n");
      fprintf(stdout,   "***************************************************\n\n");
      fprintf(stdout, "Leaks at %f Degrees\n", 180.0*atan(gamma.im/H[l].im)/M_PI);
      gammaz[l] = 0.;
      }
   else   
      gammaz[l] = 0.5*(ri[l].re*ri[l].re+ri[l].im*ri[l].im)/H[l].re;

   /* printf("gammaz[0] = %g\n",gammaz[0]); */
   /* printf("gammaz[L-1] = %g\n",gammaz[L-1]); */

   for( l = 1; l < (L - 1); l++ ){
      xlow = xpts[l]+delta;
      xhigh = xpts[l-1]-delta;
      gammaz[l] = simp_( fncn, xlow, xhigh, eps );
      /* printf("gammaz[%d] = %g \n",l,gammaz[l]); */
      }
   total = 0.0;
   
   for( l = 1; l < L-1; l++ ){
      total += gammaz[l];
      }

   /* fprintf(stdout,"Total Integrals = %15.10e \n", total); */
/*   
   xlow = xpts[L - 2];
   xhigh = xpts[0];
   printf(" Error = %15.5e \n", simp_(fncn,xlow,xhigh,eps)-total);
*/         
   total += gammaz[0] + gammaz[L-1];   
   srtotal = sqrt( total );
   for( l = 0; l < L; l++ ){
      gammaz[l] /= total;
      ri[l].re /= srtotal;
      ri[l].im /= srtotal;
      rip[l].re /= srtotal;
      rip[l].im /= srtotal;
      }

   fprintf(stdout, "\n Layer Confinement Factors \n");
   for( l = 0; l < L ; ++l)
      fprintf(stdout, "Gamma(%2d) = %15.4e %s\n", l, gammaz[l], &(comments[128*l]) );

   /*     check the integral of the fourier transform
    *
    *c      xlow =-1.d0
    *c      xhigh= 1.d0
    *c      call simp_(far,xlow,xhigh,1.d-2,total)
    *c      write(6,*)'far-field sum/2pi',total/tpi */
   return;
}/*END norm */ 
/*================================================================
 *
 *  this function evaluates the field intensity
 *
 *    field2 = phi*phi
 *
 *    where  x(i) .lt. u .lt. x(i+1), (x = xrpts)
 *
 *  input..
 *
 *    n = layers-1, the number of data points
 *    u = the abscissa at which the intensity is to be evaluated
 *
 *  if  u  is not in the same interval as the previous call, then a
 *  binary search is performed to determine the proper interval.
 **************************************************************/
double Phih2(double u)
{
  complx field;

  field = Phih(u);
  return( field.re*field.re+field.im*field.im );
}

complx Phih(double u)  
{
  int j, k, l;
  //static int TRUE = 1;
  static int i;
  double dx, hd;
  complx zc, zs, hx, field, t11, t12;
  
  if( u < xrpts[0] ){ /* Substrate -------------------------*/
    i = -1;
    l = NL - i - 2;
    dx = u - xrpts[0];
    hx = fcm(dx,hz[NL-1]);
    field = cmul(phih[NL-2],zexp(hx));
  } /*End Substrate */
  else if( u < xpts[0] ){ /* Central Layers ----------------*/
    if( i >= NL - 1 )
      i = 0;
    if( (u < xrpts[i]) || (u >= xrpts[i+1]) ) {      
      /*  binary search  */
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
    hx = fcm(dx,hz[l]);
    hd = fabs(hz[l].im*d[l]);
    if(hd < 10.0){
      t11 = ccos( hx );
      t12 = cdiv(csin(hx),hz[l]);
    }
    else{
      t11 = zexp(cmul(ZJ,hx));
      t12 = ZZERO;
    }
    zc = cmul(t11,phih[l]);
    zs = cmul(t12,phiph[l]);
    field = cadd(zc,zs);
  } /*End Central Layers*/
  else{ /* Superstrate ------------------------------------*/
    i = NL - 2;
    l = NL - i - 2;
    dx = u - xpts[0];
    hx = fcm(-dx,hz[0]);
    field = cmul(phih[0],zexp(hx));
  }

  //printf("%d %g   %g %g \n",l,u/ko,field.re,field.im); 
  return( field );
} 
/*======================================================================= */
void points(int ng, int n, double *d, double *xpts, double *xrpts)
{
  /********************************************************************
   * This routine determines the x values at the various layers assuming
   * the origin, x=0 occurs at the bottom of the layer called the grating
   * layer. The coordinate system is defined as follows:
   *
   * n-1     n-2       n-3            ng+1       ng            1       0
   * ---|-----------|--------|-...-|---------|========|-...|--------|-----
   * xpts(n-2)   xpts(n-3)      xpts(ng+1)xpts(ng)=0    xpts(1)  xpts(0)
   * xrpts(0)   xrpts(1)                                        xrpts(n-2)
   * 
   * Note that coordinates can be shifted as desired.
   * 
   *   The variables transferred to and from this routine are:
   * 
   *   ng >= 1 = grating or designated layer,
   *   n       = number of layers,
   *   d(l)    = thickness of lth layer,
   *   xpts(l) = x(l) as defined in the figure. (returned)
   *             this configuration allows the transformation matrix to
   *             be written as below.
   *   xrpts(l)= reversal of indices of xpts.
   *
   *   The field at x(n-2) is transformed to x(0) by the formula
   *
   *   Psi[x(0)] = T1.T2.T3...T(n-2) . Psi[x(n-2)]
   *
   *   where the Psi is the column vector defined below. For
   *   partial transformations to intermediate points it is possible
   *   to use groups of matrices as discussed below.  For example,
   *   the fields are transformed to the point x(ng) using the
   *   transformation:
   *
   *      Tl = T(ng+1) T(ng+2) ... T(n-2) T(n-1).
   *
   *   the transformation of field quantities using tl is
   *      Psi[x(ng)] = Tl . Psi[x(n-1)]
   *   where
   *      Psi[x] = col{phi(x),phip(x)}.
   * 
   *   The fields are transformed to the point x(ng-1) using the
   *   transformation:
   * 
   *      Tr= T(2) T(3) ...T(ng-1)
   *          
   *   using Tr gives field transformation:
   *      Psi[x(1)] = Tr . Psi[x(ng+1)]
   *******************************************************************/
  int l;

   xpts[ng] = 0;
   /*     to the right of the grating layer. */

   for( l = ng ; l > 0; l-- )
      xpts[l-1] = xpts[l] + d[l];
      
   /*     to the left to the grating layer. */

   for( l = ng + 1; l < (n - 1); l++ )
      xpts[l] = xpts[l-1] - d[l];
      

   /*     Calculation of the reversed vector.  Note that there are only
    *     n-1 points on the interior region. */

   for( l = 0; l < (n - 1); l++ )
      xrpts[l] = xpts[n - l - 2];

   for( l = 0; l < (n - 1); l++ )
     fprintf(stdout, "%2d %10.5f   %s\n", l, xpts[l]/ko, &(comments[128*l]) );

   return;
}
/*======================================================================= */
void IndexProfile(int NumberLayers, double Ko,
                  double *Xlocation, complx *Kappa)
{
  FILE *indexOut;
  int l;
  complx Index;

  indexOut = fopen("./Unloaded/IndexProfile.prn","w");
  if(indexOut == NULL){
    fprintf(stderr, "Could not open:./Unloaded/IndexProfile.prn\n");
    fprintf(stderr, "The Directory ./Unloaded must exist.\n");
    exit(0);
  }

  l = NumberLayers-1;
  Index = zsqrt(*(Kappa+l));
  fprintf(indexOut, "%10.5f  %10.5f\n",*(Xlocation+l-1)/Ko - 1.0, Index.re);
  fprintf(indexOut, "%10.5f  %10.5f\n",*(Xlocation+l-1)/Ko, Index.re);
  for(l = NumberLayers -2 ; l > 0; l --){
    Index = zsqrt(*(Kappa+l));
    fprintf(indexOut, "%10.5f  %10.5f\n",*(Xlocation+l)/Ko, Index.re);
    fprintf(indexOut, "%10.5f  %10.5f\n",*(Xlocation+l-1)/Ko, Index.re);
    }
  Index = zsqrt(*Kappa);
  fprintf(indexOut, "%10.5f  %10.5f\n",*Xlocation/Ko, Index.re);
  fprintf(indexOut, "%10.5f  %10.5f\n",*Xlocation/Ko + 1.0, Index.re);

  fclose(indexOut);
  return;
}

