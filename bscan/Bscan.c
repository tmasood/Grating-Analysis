#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <complx.h>
#include <jkbclapack.h>
#include <jkb.h>
#include <linear.h>
#include "Green.h"
#include "Dft.h"

/*----------Begin Znl GLOBAL Variables ---------------------*/
complx znl();

double *Xf;     /* Lateral field points */
double tauz;    /* 2 \tau_z = normalized tooth width */

double *xpts;   /* Array of layer boundaries */
double *xrpts;  /* Array of reverse layer boundaries */
double *gammaz; /* Gamma values for waveguide (no grating) */

complx *kappaz; /* Dielectric  constants of waveguide */
complx *hz;     /* Transvere wave numbers */
complx *chz;    /* Cosine of layers */
complx *shz;    /* Sine of layers */
complx *phi;    /* Psi at the boundaries */
complx *phip;   /* Psip at the boundaries */

/*----------Begin GLOBAL Variables -------------------------*/
int nfl;        /* Number of single-sided space harmonics */
int nflt;       /* Total space harmonics (2*nfl) */
int na;         /* Base points of first box */
int nb;         /* Base points of second box */
int nc;         /* Common side points */
int n1;         /* Points around first box (2*na+2*nc) */
int n2;         /* Points around second box (2*nb+2*nc) */
int nt;         /* Total points n1+n2 */
int nts;        /* nt*nt */
int quick;      /* quick calculation of determinant */
int nX;         /* Number of field points along x (~ 100)*/
int StoredRoots;/* Number of stored roots in calg */

double ak;      /* Grating Wavenumber 2*Pi/Lambda */
double del;     /* Normalized step along z */
double clam;    /* Normalized grating period */
double dd2;     /* del/2 */
double delc;    /* Normalized step along x */

int *IPIV;      /* Matrix Information for Lapack Routines */

int INFO;
int LWORK;
char NORM[1];
double ANORM;
double RCOND;
complx *WORK;
double *RWORK;

/* Eigenvalue Information (See Lapack) */
char JOBVL[1];
char JOBVR[1];
complx *zeign;  /* Vector of Eivenvalues */
complx *VL;     /* Matrix for left Eigenvectors */
complx *VR;     /* Matrix for right Eigenvectors */

char *Comments;

int NLayers;    /* Total number of layers */

int NGL;        /* Number of grating layers */
int GLtop;      /* Top grating layer */
int GLbot;      /* Bottom grating layer */
int *GL;        /* Pointer to grating layers */
double *DC;     /* Duty cycle of the lth grating layer*/
double *Disp;   /* Dispalcement of edge relative to bottom */
int *LKG;       /* Layer Dielectric fill in e.g. air */
int StartStack; /* First stack layer in substrate */

int *laynum;    /* Array of layer numbers */
int *islf;      /* Is the layer flat?     */
int *layreps;   /* Number of repeating layers. */
int *layrepsn;  /* Number of layers in repeating unit cell. */
int *BotLeak;   /* =1 if space harmonic leaks in substrate. */
int *TopLeak;   /* =1 if space harmonic leaks in superstrate. */

double *w;       /* Layer thicknesses in wavelengths */
double *d;       /* Normalized layer thicknesses */
double *nl;      /* Refractive index at beginning of layer */
double *nr;      /* Refractive index at end of layer */
double *kappal;  /* Dielectric Constant at begining of layer */
double *kappar;  /* Dielectric constant at end of layer */
double *bcept;   /* Dielectirc Intercept of sloping layer */
double *slp;     /* Slope of dielectric layer */
double *am;
double *tm;
double *Gamma;    /* layer confinement factors */

double *g1;       /* Greens functions */
double *h1;
double *g2;
double *h2;

complx *W;       /* exp( -j 2*Pi/N ) */
complx *WH;      /* exp( -j Pi/N) */
complx *PF;      /* Phase Factor WH^{-n(na-1)}, (0 \ldots nflt) */
complx *Vtop;    /* Ref to na/2; Propagation removed in fncn */
complx *Vtopp;
complx *Vbot;    /* Ref to na/2; Propagation removed in fncn */
complx *Vbotp;
complx *Evec;    /* Boundary Element Eigenvector as below*/
complx *VtopT;   /* Inverse transform of Vtop */
complx *VbotT;   /* Inverse transform of Vbot */
/* Counting reversed on top and left sides of boxes */
complx *q1;      /* Field points around 1st box */
complx *p1;      /* Out flux around 1st box */
complx *q2;      /* Field points around 2nd box */
complx *p2;      /* Out flux around 2nd box */
complx *Fgrat;   /* Finds in the grating nc*(na+nb) */
double *Maggrat; /* Man of Fgrat */
double *Phagrat; /* phase of Fgrat */
double *xgrat;   /* x locations in the gratins nc*(na+nb) */
double *zgrat;   /* z locations in the grating nc*(na+nb) */
double *xspl;    /* x locations in the grating nc+2 */
double *psre;    /* psi.re at the xspl points */
double *psim;    /* psi.im at the xspl points */
double *bre;     /* |---------------------------------------*/
double *cre;     /* | Real Spline coefficients of grat fld. */
double *dre;     /* |---------------------------------------*/
double *bim;     /* |--------------------------------------------*/
double *cim;     /* | Imaginary Spline coefficients of grat fld. */
double *dim;     /* |--------------------------------------------*/
/*-----------------------------------------------------------------
 * Unprimes are fields and primes are flux
 * -------------------------------------------------
 *|         2           |             5             |
 *|                     |nc                         |
 *|                    1|6'                        4|7'
 *|                     |                           |
 *|       0   na        |         3  nb             |
 * -------------------------------------------------
 *----------------------------------------------------------------*/
complx *z;   /* Boundary element storage matrix */
complx *h;   /* Transverse wavenumbers index by n_+l*nflt */
complx *c;   /* cos(h*d) */
complx *s;   /* sin(h*d) */
complx *t;   /* tan(h*d) */
complx *r;   /* czeta tan(h*d) - j */
complx *cm;  /* psip/psi at bottom of grating */
complx *cn;  /* psip/psi at bottom of grating (needs a negative sign) */
complx *tcm; /* Transform of cm */
complx *tcn; /* Transform of cn */
complx *zm;  /* tcm with propagation terms? Nevertheless p = zm*u */
complx *zn;  /* tcn with propagation terms? Nevertheless p = zn*u */
complx *zgam;/* exp( -gamma z ) */
complx *An;  /* psi at layer boundaries (Amplitude of Cosine)  Nflt*Layers */
complx *Bn;  /* pipp at layer boundaries (Amplitude of Sin/h)   Nflt*Layers */

complx *phin;/* Working vector for phase modified fields, n = 0, \ldot, nflt. */
complx *Phil;/* Working vector for z field values, l = 0, \ldot, nflt. */

int iflag, nrts;

double cl[10], al[10], be[10];

FILE *fin;
FILE *fMagFields;
FILE *fPhaFields;
FILE *fLoadedField;
FILE *fOmegaAlpha;
FILE *fOmegaBeta;
/*-----------BEGIN GLOBAL FUNCTIONS---------------------------------*/
void ReadError(int Line, char *outstring );

int readdata(double *lamum, double *clamum, double *dc,
             int *nfl, int *nc, 
	     double *zguessr, double *zguessi,
             double *varend, double *dvar, double *pecmax, 
             int *itmin, int *itmax,
	     int *NL, int *NGL, int *startstack );

void ReadLayers( int *laynum, int *islf, int *layreps, int *layrepsn,
		 double *nl, double *nr, double *w,
		 char *Comments, int N, int Line);

int ReadGratLayers( int ,int , int *, double *, double *, 
		    int *, int );

void calgs(int store,  int n,  int *ihaf, 
	   int *idub, double pecmax, complx *z1, complx *z2, complx *z3, 
	   double *c, double *dc);

void InitializeData(int N, double lambda, double *ko, double *nl, double *nr,
		    double *kl, double *kr, double *w, double *d, 
		    double *bcept, double *slp, double *am, double *tm);

complx fncn(complx );

void GetWaveguideParameters(int , int *, int *, 
			    double *, double *, complx *);

void floquet( int Nt, int Nf, int Ntop, int Nbot, int NL,
	      complx *W, complx *WH, complx *Vt, complx *Vb, 
	      complx *h, complx *c, complx *s, complx *t,
	      complx *A, complx *B );

double Norm(int nlayers, int nfl, int nflt, int ntop, int nbot,
	    char *comments,
	    int *TopLeak, int *BotLeak,
	    double *dum, double *Xpts, double *Xrpts, double *zgrat, 
	    double tauz, double *Gam,
	    complx *ht, complx *ct, complx *st, complx *tt,
	    complx *ps, complx *psp, complx *W, complx *PF,
	    int nspl,
	    double *xspl_, double *psre_, double *psim_,
	    double *bre_, double *cre_, double *dre_,
	    double *bim_, double *cim_, double *dim_);

void ModeSketch(int nitter, int nX, double *Xf, double ko,
		complx *phin, complx *Phil);

void GratingSketch(int nitter, int nflt, int na, int nb, int nc, 
		   double delc, double del, double fko,
		   int GLbot, int GLtop, double *xpts,
		   double *xgrat, double *zgrat, complx *Fgrat,
		   complx g, double ng, double nfill,
		   complx *Vtop, complx *Vbot,
		   complx *p1, complx *q1, complx *p2, complx *q2 );

complx Psi(double );


/*------------END GLOBAL FUNCTIONS---------------------------------*/

int main( int argc, char *argv[] )
{
  int store;
  int i, idub, ihaf, ingap, il, 
    itmax, itmin, itt, l,
    nitter;
  int Nxgrating;
  int Nzgrating;
  int Nzstart, Nzend;
  double Vtopmag, Vtoppha, Vbotmag, Vbotpha;
  double Fnorm;

  int *ix1;  /* x center points of A panels */
  int *iz1;  /* z center points of A panels */
  int *ix2;  /* x center points of B panels */
  int *iz2;  /* z center points of B panels */
  int *ir1;  /* Related to |r1-r2| in A */
  int *ir2;  /* Related to |r1-r2| in B */
  double lamum, alcm, alpha, beta, betan, dc, dccal, delcum, 
    delta, delum, dvar, fko, pecmax, pectemp, varend, 
    zgi, zgr;
  complx g, z1, z2, z3, zroot, 
    zw;
  complx Zg[3];

  /* EQUIVALENCE translations */
  double _e0[1];
  double *const clamum = (double*)_e0;
  double *const var = (double*)_e0;
  /* end of EQUIVALENCE translations */

  int line;
  int cLine, CalFields = 0, UseGuess = 0;

  StoredRoots = 0;

  while ( --argc > 0 && ( *++argv)[0] == '-' )
    while ( cLine = *++argv[0] )
      switch (cLine) {
      case 'F':
	CalFields = 1;
	break;
      case 'G':
	UseGuess = 1;
	break;
      default:
	printf("Illegal Option: %c\n", cLine);
	argc = 0;
	exit( EXIT_FAILURE );
	break;
      }

  if(argc != 1 ) {
    fprintf(stderr,"\nUsage:\tBscan -F -G <input file> \n\n");
    fprintf(stderr,"\t-F option to calculate Fields.\n");
    fprintf(stderr,"\t-G option to used Guessed Root.\n");
    fprintf(stderr,"\n");
    exit( EXIT_FAILURE );
  }
  else {
    fin = fopen( *argv, "r");
    if(fin == NULL){
      fprintf(stderr, "Could not open '%s'\n", *argv );
      exit( EXIT_FAILURE );
    }
  }
  /*-------------------------------------------------------------*/
  fOmegaAlpha = fopen("OMB/alpha", "w");
  if(fOmegaAlpha == NULL){
    fprintf(stderr, "Could not open 'OMB/alpha'\n");
    exit( EXIT_FAILURE );
  }

  fOmegaBeta = fopen("OMB/beta", "w");
  if(fOmegaBeta == NULL){
    fprintf(stderr, "Could not open 'OMB/beta'\n");
    exit( EXIT_FAILURE );
  }

  line = readdata( &lamum, clamum, &dc, &nfl, &nc, &zgr, &zgi, 
		   &varend, &dvar, &pecmax, &itmin, &itmax, 
		   &NLayers, &NGL, &StartStack );

  /*-------Start Allocation of Layer Arrays ---------------------*/
  Comments = (char *)calloc(NLayers*128,sizeof(char));

  laynum = (int *)calloc(NLayers,sizeof(int));
  islf = (int *)calloc(NLayers,sizeof(int));
  layreps = (int *)calloc(NLayers,sizeof(int));
  layrepsn = (int *)calloc(NLayers,sizeof(int));

  w = (double *)calloc(NLayers,sizeof(double));
  d = (double *)calloc(NLayers,sizeof(double));
  nl = (double *)calloc(NLayers,sizeof(double));
  nr = (double *)calloc(NLayers,sizeof(double));
  kappal = (double *)calloc(NLayers,sizeof(double));
  kappar = (double *)calloc(NLayers,sizeof(double));
  bcept = (double *)calloc(NLayers,sizeof(double));
  slp = (double *)calloc(NLayers,sizeof(double));
  am = (double *)calloc(NLayers,sizeof(double));
  tm = (double *)calloc(NLayers,sizeof(double));
  Gamma = (double *)malloc(NLayers*sizeof(double));

  GL = (int *)calloc(NLayers,sizeof(int));
  Disp = (double *)calloc(NLayers,sizeof(double)); /* For Furture */
  DC = (double *)calloc(NLayers,sizeof(double));
  LKG = (int *)calloc(NLayers,sizeof(int));

  xpts   = (double *)calloc(NLayers,sizeof(double));
  xrpts  = (double *)calloc(NLayers,sizeof(double));
  Xf = (double *)calloc(100,sizeof(double)); /* Also Nfld in znl */
  gammaz = (double *)calloc(NLayers,sizeof(double));

  kappaz = (complx *)calloc(NLayers,sizeof(complx));
  hz = (complx *)calloc(NLayers,sizeof(complx));
  chz = (complx *)calloc(NLayers,sizeof(complx));
  shz = (complx *)calloc(NLayers,sizeof(complx));
  phi = (complx *)calloc(NLayers,sizeof(complx));
  phip = (complx *)calloc(NLayers,sizeof(complx));

  /*-------End Allocation of Layer Arrays ---------------------*/
  nflt = 2*nfl;
  na = dc*nflt ;
  if(!(na & 1)){
    fprintf(stderr, "Boundary points under the tooth must be odd!\n");
    fprintf(stderr, "... na = %d\n",na);
    exit(EXIT_FAILURE);
  }
  nb = nflt - na;
  n1 = 2*(na + nc);
  n2 = 2*(nb + nc);
  nt = n1 + n2;
  nts = nt*nt;
  nX = 100;
  /*------Start Allocation of Matrix and Vector Arrays --------*/
  ix1 = (int *)malloc(n1*sizeof(int));
  iz1 = (int *)malloc(n1*sizeof(int));
  ix2 = (int *)malloc(n2*sizeof(int));
  iz2 = (int *)malloc(n2*sizeof(int));
  ir1 = (int *)malloc(n1*n1*sizeof(int));
  ir2 = (int *)malloc(n2*n2*sizeof(int));

  TopLeak = (int *)malloc(nflt*sizeof(int));
  BotLeak = (int *)malloc(nflt*sizeof(int));

  W =  (complx *)malloc(nflt*sizeof(complx));
  WH = (complx *)malloc(2*nflt*sizeof(complx));
  PF = (complx *)malloc(nflt*sizeof(complx));
  Vtop = (complx *)calloc(nflt,sizeof(complx));
  Vtopp = (complx *)calloc(nflt,sizeof(complx));
  Vbot = (complx *)calloc(nflt,sizeof(complx));
  Vbotp = (complx *)calloc(nflt,sizeof(complx));
  VtopT = (complx *)malloc(nflt*sizeof(complx));
  VbotT = (complx *)malloc(nflt*sizeof(complx));
  phin = (complx *)malloc(nflt*sizeof(complx));
  Phil = (complx *)malloc(nflt*sizeof(complx));
  Fgrat = (complx *)malloc(nc*(na+nb)*sizeof(complx));
  Maggrat = (double *)malloc(nc*(na+nb)*sizeof(double));
  Phagrat = (double *)malloc(nc*(na+nb)*sizeof(double));
  xgrat = (double *)malloc(nc*(na+nb)*sizeof(double));
  zgrat = (double *)malloc(nc*(na+nb)*sizeof(double));
  xspl = (double *)malloc((nc+2)*sizeof(double));
  psre = (double *)malloc((nc+2)*sizeof(double));
  psim = (double *)malloc((nc+2)*sizeof(double));
  bre = (double *)malloc((nc+2)*sizeof(double));
  cre = (double *)malloc((nc+2)*sizeof(double));
  dre = (double *)malloc((nc+2)*sizeof(double));
  bim = (double *)malloc((nc+2)*sizeof(double));
  cim = (double *)malloc((nc+2)*sizeof(double));
  dim = (double *)malloc((nc+2)*sizeof(double));

  LWORK = nts;
  IPIV = (int *)calloc(nt,sizeof(int));
  WORK =  (complx *)malloc(LWORK*sizeof(complx)); /* 4 instead of 2 */
  RWORK = (double *)malloc(LWORK*sizeof(double)); /* 4 instead of 2 */
  NORM[0] = '1';
  JOBVL[0] = 'V';
  JOBVR[0] = 'V';

  cm = (complx *)malloc(nflt*sizeof(complx));
  cn = (complx *)malloc(nflt*sizeof(complx));
  tcm = (complx *)malloc(nflt*sizeof(complx));
  tcn = (complx *)malloc(nflt*sizeof(complx));
  zm = (complx *)malloc(nflt*nflt*sizeof(complx));
  zn = (complx *)malloc(nflt*nflt*sizeof(complx));
  zgam = (complx *)malloc(nflt*sizeof(complx));

  h = (complx *)malloc(nflt*NLayers*sizeof(complx ));
  c = (complx *)malloc(nflt*NLayers*sizeof(complx ));
  s = (complx *)malloc(nflt*NLayers*sizeof(complx ));
  t = (complx *)malloc(nflt*NLayers*sizeof(complx ));
  r = (complx *)malloc(nflt*NLayers*sizeof(complx ));
  An = (complx *)malloc(nflt*NLayers*sizeof(complx));
  Bn = (complx *)malloc(nflt*NLayers*sizeof(complx));
  z = (complx *)calloc(nts,sizeof(complx));
  zeign = (complx *)malloc(nt*sizeof(complx));
  Evec = (complx *)malloc(nt*sizeof(complx));
  VR = (complx *)malloc(nts*sizeof(complx));
  VL = (complx *)malloc(nts*sizeof(complx));
  q1 = (complx *)malloc(n1*sizeof(complx ));
  p1 = (complx *)malloc(n1*sizeof(complx ));
  q2 = (complx *)malloc(n2*sizeof(complx ));
  p2 = (complx *)malloc(n2*sizeof(complx ));
  g1 = (double *)calloc(n1*n1,sizeof(double));
  h1 = (double *)calloc(n1*n1,sizeof(double));
  g2 = (double *)calloc(n2*n2,sizeof(double));
  h2 = (double *)calloc(n2*n2,sizeof(double));
  /*--------End Allocation of Matrix and Vector Arrays --------*/

  line = ReadGratLayers( NLayers, NGL, GL, Disp, DC, LKG, line);

  ReadLayers( laynum, islf, layreps, layrepsn, nl, nr, w, 
	      Comments, NLayers, line);

  fclose(fin);
  /*------Calculate Array Parameters --------------------------*/
  mfft( nflt, W );
  mfftH( nflt, WH );
  PhaseFactor(nflt, nfl, na, WH, PF);

  for(i = 0; i < NLayers; i++){
    GLtop = i;
    if(GL[i])
      break;
  }
  for(i = NLayers-2; i > 0; i--){
    GLbot = i;
    if(GL[i])
      break;
  }
  pectemp = pecmax;
  ingap = 0;
  nrts = 0;
  nitter = 0;
  idub = 0;
  ihaf = 0;
  store = 0;
  InitializeData(NLayers, lamum, &fko, nl, nr, kappal, kappar, w, d, 
		 bcept, slp, am, tm);

  printf(" lambda =  %g \n", lamum );
  printf("\n layer   index      thickness\n");
  printf("------------------------------------\n");
  for( i = 0; i < NLayers ; i++ )
    printf("%3d %10.6f %10.6f %10.6f  %s\n", laynum[i], nl[i], nr[i], w[i], 
	   &(Comments[128*i]) );
  printf("------------------------------------\n");

  /* Initialize Kappaz */
  GetWaveguideParameters(NLayers, GL, LKG, DC, kappal, kappaz);

  z3 = ftoc( zgr, zgi );

  /* Return the unloaded root. */
  z3 = znl(Comments, itmax, NLayers, GL, DC, fko, d, xrpts, xpts, Xf, gammaz,
	   kappaz, hz, chz, shz, phi, phip, z3);

  if( UseGuess )
    z3 = ftoc( zgr, zgi );

  /*-----------------------------------------------------------
   * Do a scan about a given Bragg Wavelength:
   * clamum = (2pi/ko / n_eff) * 0.5 for 1st Bragg
   * clamum = (2pi/ko / n_eff) * 1.0 for 2nd Bragg
   *----------------------------------------------------------*/

  /* for First Bragg */
  //  *clamum = 0.5*lamum/z3.im;
  //varend = 1.1*(*clamum); /* 1 percent above */
  //*var = 0.8*(*clamum);
  //  dvar = (varend - *var)/10;
  
  //varend = *var;
  //*clamum = 0.4;
  //varend = 0.3;
  //dvar = -0.01;

  *var = *var - dvar; /* Start with initial calculated value of var */

  //while( *var > varend ) {
  while( nitter < 1 ) {
    calgs( store, StoredRoots, &ihaf, &idub, pecmax, &z1, &z2, &z3, 
	   var, &dvar );

    dccal = (double )nflt;
    dccal = (double )na/dccal;
    printf( "***************************************************************\n");
    printf( "***   Itteration %d \n", nitter );
    printf( "***************************************************************\n");
    printf( "Panel Sizes %d %d %d \n", na, nb, nc );
    printf( "Grating Heigt %g \n", w[1] );
    printf( "Calculated Duty Cycle %f \n", dccal );

    /* Itteration Loop starts here */
    delum = *clamum/(double )nflt;
    printf( "Tooth Size (A region)  %g \n", na*delum );
    printf( "Fill Size  (B region)  %g \n", nb*delum );
    printf( "Grating Period = %g\n", nflt*delum);
    printf( "Boundary Matrix N*N, N = %d\n\n", nt);

    delcum = w[1]/(double )nc;
    del = fko*delum;
    dd2 = 0.5*del;
    tauz = na*dd2;
    delc = fko*delcum;
    clam = nflt*del;
    ak = lamum/(*clamum);

    boupts( na, nb, nc, n1, n2, iz1, iz2, ix1, 
	    ix2, ir1, ir2 );
    /*
      printf("iz1\n");
      iout_(5,5,iz1);
      printf("ix1\n");
      iout_(5,5,ix1);
      printf("iz2\n");
      iout_(5,5,iz2);
      printf("ix2\n");
      iout_(5,5,ix2);

      printf("ir1\n");
      iout_(n1,n1,ir1);
    */
    green( na, nc, n1, iz1, ix1, ir1, del, delc, nl[1], g1, h1 );
    dodiag( na, nc, n1, nl[1], del, delc, g1, h1 );
    dorev( na, nc, n1, g1, h1);
    printf("Completed g1 and h1...\n");
    green( nb, nc, n2, iz2, ix2, ir2, del, delc, nl[0], g2, h2 );
    dodiag( nb, nc, n2, nl[0], del, delc, g2, h2 );
    dorev( nb, nc, n2, g2, h2 );
    printf("Completed g2 and h2...\n");
    /*
      printf("g2\n");
      rout_(n2, n2, g2);
      printf("h2\n");
      rout_(n2, n2, h2);
    */
    //  z3 = ftoc( zgr, zgi );
    //  zw = fncn(z3);
    Zg[0] = z1;
    Zg[1] = z2;
    Zg[2] = z3;

    iflag = 0;
    quick = 1;
    printf(" Minimizing the determinant...\n");
    zroot = croot( fncn, Zg, &zw, &itt, itmax, &store );
    /* Decide on a Root */
    if (store){
      itmax = itmin; /* Root search should be quick */
      printf(" Minimizing the smallest eigenvalue...\n");
      quick = 0;
      for(i = 0; i < 3 ; i++)
	Zg[i].re = fabs(Zg[i].re);
      zroot = croot( fncn, Zg, &zw, &itt, itmax, &store );
      ++StoredRoots;
      fprintf(fOmegaAlpha,"%f\t%f\n",clam*zroot.re, clam);
      fprintf(fOmegaBeta,"%f\t%f\n",clam*zroot.im/M_2PI, clam);
      z3 = zroot;
      g = zroot;
      alpha = g.re;
      beta = g.im;
      betan = fko*beta*(*clamum);
      alcm = fko*alpha;
      printf( " Alpha = %e \n", alcm );
      printf( " Beta*Clambda/2Pi = %f\n", betan/M_2PI );
      if(CalFields){
	printf(" Calculate the Eigenvectors...\n");
	/* Calculate the Eigenvector */
	quick = -1;
	zw = fncn(zroot);

	/* Calculate the fields (not normalized) at center of tooth */
	Nzstart = na >> 1;
	Nzend = Nzstart + 1;
	xspl[0] = 0.0;
	xspl[nc+1] = delc*nc;
	for( l = Nzstart; l < Nzend ; l++)
	  zgrat[l] = 0.0;
	for( i = 0; i < nc ; i++){
	  xgrat[i] = (2*i+1)*delc*0.5;
	  xspl[i+1] = xgrat[i];
	}
	Nxgrating = nc;
	Nzgrating = Nzend-Nzstart;
	GratingFields(na, nc, nb, del, delc,p1, q1, p2, q2,
		      Nxgrating, Nzgrating, xgrat, zgrat, 
		      g, Fgrat, nl[1], nl[0]);
	psre[0] = Vbot[Nzstart].re;
	psim[0] = Vbot[Nzstart].im;
	psre[nc+1] = Vtop[Nzstart].re;
	psim[nc+1] = Vtop[Nzstart].im;
	for(i = 0; i < nc ; i++){
	  psre[i+1] = Fgrat[i].re;
	  psim[i+1] = Fgrat[i].im;
	}
	/* spline(nc+2, xspl, psre, bre, cre, dre);
	spline(nc+2, xspl, psim, bim, cim, dim); */

	floquet(nflt, nfl, GLtop, GLbot, NLayers, W, WH, Vtop, Vbot,
		h, c, s, t, An, Bn);
	//  exit( 0 );

	Fnorm =  Norm( NLayers, nfl, nflt, GLtop, GLbot, Comments,
		       TopLeak, BotLeak, d, xpts, xrpts, zgrat, tauz,
		       Gamma, h, c, s, t, An, Bn, W, PF, nc+2,
		       xspl, psre, psim, bre, cre, dre, bim, cim, dim);

	// Normalize the fields in the grating 
	for( i = 0; i < nflt; i++) {
	  Vtop[i] = fcm( Fnorm,Vtop[i] );
	  Vbot[i] = fcm( Fnorm,Vbot[i] );
	}

	for( i = 0; i < n1; i++) {
	  p1[i] = fcm( Fnorm, p1[i] );
	  q1[i] = fcm( Fnorm, q1[i] );
	  p2[i] = fcm( Fnorm, p2[i] );
	  q2[i] = fcm( Fnorm, q2[i] );
	}
	// End Grating field normalization

	GratingSketch(nitter, nflt, na, nb, nc, delc, del, fko,
		      GLbot, GLtop, xpts, xgrat, zgrat, Fgrat,
		      g, nl[GLbot], nl[GLtop-1], Vtop, Vbot,
		      p1, q1, p2, q2);

	ModeSketch(nitter, nX, Xf, fko, phin, Phil );
      }
    }
    else {
      ++ihaf;
      *var -= dvar;
      dvar *= 0.5;
      printf(" No Convergence: Interval Halved.\n");
      idub = 0;
      if(ihaf > 4 ) {
	printf("4 Consecutive interval halves.\n");
	exit( EXIT_FAILURE );
      }
    }
    ++nitter;
  }
  return( EXIT_SUCCESS );
}
/*------------------------------------------------------------------------*/
void GetWaveguideParameters(int NL, int *GL, int *Kf, double *DC, 
			    double *K, complx *Kz)
{
  int l;
  
  for( l = 0; l < NL; l++){
    if(GL[l]){
      Kz[l].re = K[l]*DC[l] + K[Kf[l]]*(1.0-DC[l]);
      Kz[l].im = 0.0;
    }
    else{
      Kz[l] = ftoc(K[l],0.0);
    }
  }
  return;
}
/*------------------------------------------------------------------------*/
void ReadError(int Line, char *outstring )
{
  fprintf(stderr, "Read error on Line %d\n", Line);
  fprintf(stderr, "Input String was: %s\n", outstring );
  exit( EXIT_FAILURE );
  return;
}
/*----------------------------------------------------------------------- */
int readdata(double *lamum, double *clamum, double *dc,
             int *nfl, int *nc, 
	     double *zguessr, double *zguessi,
             double *varend, double *dvar, double *pecmax, 
             int *itmin, int *itmax,
	     int *NL, int *NGL, int *startstack )
{
  char buffer[80];
  int nbuf = 80;
  int line;

  //  *itmin = 10;
  // *itmax = 50;

  line = 0;

  while(fgets(buffer,nbuf,fin) != NULL ){
    line++;
    if(*buffer != '#'){
      if( sscanf(buffer, "%lf", lamum ) < 1 )
	ReadError(line, buffer);
      break;
    }
  }

  while(fgets(buffer,nbuf,fin) != NULL ){
    line++;
    if(*buffer != '#'){
      if( sscanf(buffer, "%lf %lf",  clamum, dc ) < 2 )
	ReadError(line, buffer);
      break;
    }
  }

  while( fgets(buffer,nbuf,fin) != NULL ){
    line++;
    if(*buffer != '#'){
      if( sscanf(buffer, "%d %d", nfl, nc ) < 2 ) 
	ReadError(line, buffer); 
      break;
    }
  }

  while( fgets(buffer,nbuf,fin) != NULL ){
    line++;
    if(buffer[0] != '#'){
      if( sscanf( buffer, "%lf %lf", zguessr, zguessi ) != 2 )
	ReadError(line, buffer);
      break;
    }
  } 

  while( fgets(buffer,nbuf,fin) != NULL ){
    line++;
    if(buffer[0] != '#'){
      if( sscanf( buffer, "%d %d", itmin, itmax ) != 2 )
	ReadError(line, buffer);
      break;
    }
  } 

  while( fgets(buffer,nbuf,fin) != NULL ){
    line++;
    if(buffer[0] != '#'){
      if( sscanf( buffer, "%lf %lf %lf", varend, dvar, pecmax ) < 3 )
	ReadError(line, buffer);
      break;
    }
  }

  while( fgets(buffer,nbuf,fin) != NULL ){
    line++;
    if(buffer[0] != '#'){
      if( sscanf( buffer, "%d %d %d", NL, NGL, startstack ) < 3 )
	ReadError(line, buffer);
      break;
    }
  } 

  return( line );
}
/*---------------------------------------------------------------------------*/
void ReadLayers( int *laynum, int *islf, int *layreps, int *layrepsn,
		 double *nl, double *nr, double *w, 
		 char *Comments, int N, int Line)
     /**************************************************************
      * The data file must be arranged as follows
      * (# represents a comment on the line.)
      *
      * layer_number   index  	width  	  comments
      * .
      * .
      * .
      **************************************************************/
{
  char buffer[80];
  int nbuf = 80;
  int l;
  double dummy;

  l = 0;
  while( l < N ){
    while( fgets(buffer,nbuf,fin) != NULL ){
      Line++;
      if(*buffer != '#'){
        if(sscanf(buffer, "%d %lf %lf %lf %d %lf %d %d %s", laynum+l, 
		  nl+l, nr+l, &dummy, islf+l, w+l, layreps+1, layrepsn+l,
		  Comments+l*128 ) < 9 )
          ReadError(Line, buffer);
        break;
      }
    }
    l++;
  }

  return;
}
/*---------------------------------------------------------------------*/
int ReadGratLayers( int NL, int L, int *GL, 
		    double *Displacement, double *DutyCycle, 
		    int *LKG, int Line)
{
  char buffer[80];
  int nbuf;
  int l;
  int lay;
  int K;
  double dc;
  double disp;

  for (l = 0; l < NL; l++){
    DutyCycle[l] = 1.0;
    LKG[l] = l;
  }

  l = 0;
  while( l < L ){
    while( fgets(buffer,nbuf,fin) != NULL ){
      Line++;
      if(*buffer != '#'){
        if(sscanf(buffer, "%d %lf %lf %d", &lay, &disp, &dc, &K) < 4)
	  ReadError(Line, buffer);
	*(GL+lay) = lay;
	*(Displacement+lay) = disp;
	*(DutyCycle+lay) = dc;
	*(LKG+lay) = K;
        break;
      }
    }
    l++;
  }
  return( Line );
}
/*---------------------------------------------------------------------*/
void InitializeData(int N, double lambda, double *ko, double *nl, double *nr,
		    double *kl, double *kr, double *w, double *d, 
		    double *bcept, double *slp, double *am, double *tm)
{
  int l;
  *ko = 2.0*M_PI/lambda;

  for( l = 0 ; l < N ; l++ ){
    kl[l] = nl[l]*nl[l];
    kr[l] = nr[l]*nr[l];
    d[l] = *ko*w[l];
    if(nl[l] != nr[l]){
      bcept[l] = kl[l];
      slp[l] = (kr[l] - kl[l])/d[l];
      am[l] = cbrt(slp[l]*slp[l]);
      tm[l] = -copysign( cbrt(fabs(slp[l])),slp[l] );
    }
  }
  return;
}
/*---------------------------------------------------------------------*/
complx TopTrans(int n_)
/* Transformation across the bottom layers for n_ space harmonic. */
{
  int l, nl;
  complx vtop, vtopp, r1, r2;
  complx psipOpsi;

  vtop.re = 1.0;
  vtop.im = 0.0;
  vtopp = cneg(h[n_]);
  
  for(l = 1; l < GLtop; l++){
    nl = n_ + l*nflt;
    r1 = csub(vtop,tdiv(vtopp,t[nl],h[nl])); 
    r2 = cadd(vtopp,tmul(vtop,t[nl],h[nl]));
    vtop = r1;
    vtopp= r2;
  }
  psipOpsi = cdiv(vtopp,vtop);
  return( psipOpsi);
}
/*---------------------------------------------------------------------*/
complx BotTrans(int n_)
/* Transformation across the bottom layers for n_ space harmonic. */
{
  int l, nl;
  complx vbot, vbotp, b1, b2;
  complx psipOpsi;

  l = NLayers-1;
  vbot.re = 1.0;
  vbot.im = 0.0;
  vbotp = h[n_+l*nflt];
  for(l = NLayers-2; l > GLbot; l--) {
    nl = n_ + l*nflt;
    b1 = cadd(vbot,tdiv(vbotp,t[nl],h[nl]));
    b2 = csub(vbotp,tmul(vbot,t[nl],h[nl]));
    vbot = b1;
    vbotp = b2;
  }
  psipOpsi = cdiv(vbotp,vbot);
  return( psipOpsi);
}
/*---------------------------------------------------------------------*/
complx fncn(complx g)
{
  int i, ii, ij, ijg, ijh, ijz, il, j, l, lj, n, n_, nl;
  int llp, lp, L, LA;
  int perms;
  int IMIN;
  int N0, N1, N2, N3, N4, N5, N6, N7, N8;
  int icenter;
  int OptWork;

  double dnflt;
  double rval;
  double RMIN;
  double Vbotmag, Vbotpha, Vtopmag, Vtoppha;

  complx det, 
    ex, gn;

  complx gn2;
  complx ZMIN;
  complx vcenter;

  //  void zgedi(), zgefa();

  //  static complx J = {0.e0,1.e0};
  static complx Z = {0.0,0.0};

  dnflt = (double )nflt;
  ++nrts;
  for( n = -nfl; n < nfl ; n++ ){
    gn.re = g.re;
    gn.im = g.im+n*ak;
    n_ = nfl+n;
    gn2 = csq( gn );
    l = 0;
    nl = n_ + l*nflt;
    h[nl] = zsqrt( cneg(cfadd(gn2,kappal[l])) );
    for(l = 1; l < NLayers-1; l++){
      if(!GL[l]){
	nl = n_ + l*nflt;
	h[nl] = zsqrt( cfadd(gn2,kappal[l]) );
	c[nl] = ccos(cfm(h[nl],d[l]));
	s[nl] = csin(cfm(h[nl],d[l]));
	t[nl] = ztan(cfm(h[nl],d[l]));
	r[nl] = czeta(cfm(h[nl],d[l]));
      }
    }
    l = NLayers-1;
    nl = n_ + l*nflt;
    h[nl] =  zsqrt( cneg(cfadd(gn2,kappal[l])) ); 

    // Starting New Code
    cn[n_] = cdd(TopTrans(n_),dnflt);
    cm[n_] = cdd(BotTrans(n_),-dnflt); /* n is in neg x dir */
  }
  //  printf("cm\n");
  //  cout_(nflt,1,cm);
  /* Calculate ZM and ZN Matrices */
  dft(nflt, nfl, W, cm, tcm);
  dft(nflt, nfl, W, cn, tcn);

  /* zgam is used only for differencing 
   * i.e., zgam_{mn} = exp(-\gamma (z_m-z_n)) */
  zgam[0] = ftoc(1.0,0.0);
  for(l = 1; l < nflt; l++){
    zgam[l] = zexp(cfm(g, -l*del));
  }

  for( l = 0; l < nflt ; l++ ){
    for(lp = 0; lp < nflt; lp++ ){
      llp = l+nflt*lp;
      L = l-lp;
      if(L < 0 ){
	LA = abs(L);
       	zm[llp] = cdiv(tcm[L+nflt],zgam[LA]);
       	zn[llp] = cdiv(tcn[L+nflt],zgam[LA]);
      }
      else{
       	zm[llp] = cmul(tcm[L],zgam[L]);
       	zn[llp] = cmul(tcn[L],zgam[L]);
      }
    }
  }

  /*
    printf("ZM\n");
    cout_(nflt,nflt,zm);
    printf("ZN\n");
    cout_(nflt,nflt,zn);
  */

  /*     BEGIN LOADING THE Z MATRIX
   ****ELEMENT 11
   * */
  for( j = 0; j < na; j++ ){
    for( i = 0; i < n1; i++ ){
      ij = i + j*nt;
      ijh = i + j*n1;
      z[ij] = ftoc(h1[ijh],0.);
      for( l = 0; l < na; l++ ){
	il = i + l*n1;
	lj = l + j*nflt;
	z[ij] = csub(z[ij],fcm(g1[il],zm[lj]));
      }
    }
  }
  /****ELEMENT 21
   * */
  for( j = 0; j < na; j++ ){
    for( i = n1; i < nt; i++ ){
      ij = i + j*nt;
      z[ij] = Z;
      for( l = 0; l < nb; l++ ){
	il = i - n1 + l*n2;
	lj = l + na + j*nflt;
	z[ij] = csub(z[ij],fcm(g2[il],zm[lj]));
      }
    }
  }

  /****ELEMENT 12
   * */
  for( j = na; j < na + nc; j++ ){
    for( i = 0; i < n1; i++ ){
      ijz = i + j*nt;
      ijh = i + j*n1;
      z[ijz] = ftoc(h1[ijh],0.);
    }
  }

  /****ELEMENT 22
   * */
  for( j = na; j < na + nc; j++ ){
    for( i = n1; i < nt; i++ ){
      ijz = i + j*nt;
      ijh = i - n1 + (j + 2*nb + nc - na)*n2;
      z[ijz] = ftoc(h2[ijh],0.);
    }
  }

  /****ELEMENT 13
   * */
  for( j = na + nc; j < 2*na + nc; j++ ){
    for( i = 0; i < n1; i++ ){
      ij = i + j*nt;
      ijh = i + j*n1;
      z[ij] = ftoc(h1[ijh],0.);
      for( l = 0; l < na; l++ ){
	il = i + (l + na + nc)*n1;
	lj = l + (j - na - nc)*nflt;
	z[ij] = csub(z[ij],fcm(g1[il],zn[lj]));
      }
    }
  }

  /****ELEMENT 23 */
  for( j = na + nc; j < 2*na + nc; j++ ){
    for( i = n1; i < nt; i++ ){
      ij = i + j*nt;
      z[ij] = Z;
      for( l = 0; l < nb; l++ ){
	il = i - n1 + (l + nb + nc)*n2;
	lj = l + na + (j - na - nc)*nflt;
	z[ij] = csub(z[ij],fcm(g2[il],zn[lj]));
      }
    }
  }
  
  /****ELEMENT 14
   * */
  for( j = 2*na + nc; j < 2*na + nc + nb; j++ ){
    for( i = 0; i < n1; i++ ){
      ij = i + j*nt;
      z[ij] = Z;
      for( l = 0; l < na; l++ ){
	il = i + l*n1;
	lj = l + (j - na - nc)*nflt;
	z[ij] = csub(z[ij],fcm(g1[il],zm[lj]));
      }
    }
  }

  /****ELEMENT 24
   * */
  for( j = 2*na + nc; j < 2*na + nc + nb; j++ ){
    for( i = n1; i < nt; i++ ){
      ij = i + j*nt;
      ijh = i - n1 + (j - 2*na - nc)*n2;
      z[ij] = ftoc(h2[ijh],0.);
      for( l = 0; l < nb; l++ ){
	il = i - n1 + l*n2;
	lj = l + na + (j - na - nc)*nflt;
	z[ij] = csub(z[ij],fcm(g2[il],zm[lj]));
      }
    }
  }

    /****ELEMENT 15
     * */
    ex = zexp( cfm(g,clam) );
    for( j = 2*na + nb + nc; j < n1 + nb; j++ ){
      for( i = 0; i < n1; i++ ){
	ijz = i + j*nt;
	ijh = i + (j - nb)*n1;
	z[ijz] = cfm(ex,h1[ijh]);
      }
    }

    /****ELEMENT 25
     * */
    for( j = 2*na + nb + nc; j < 2*na + nb + 2*nc; j++ ){
      for( i = n1; i < nt; i++ ){
	ijz = i + j*nt;
	ijh = i - n1 + (j - 2*na - nc)*n2;
	z[ijz] = ftoc(h2[ijh],0.);
      }
    }

  /****ELEMENT 16
   * */
  for( j = n1 + nb; j < n1 + 2*nb; j++ ){
    for( i = 0; i < n1; i++ ){
      ij = i + j*nt;
      z[ij] = Z;
      for( l = 0; l < na; l++ ){
	il = i + (l + na + nc)*n1;
	lj = l + (j + na - n1 - nb)*nflt;
	z[ij] = csub(z[ij],fcm(g1[il],zn[lj]));
      }
    }
  }

  /****ELEMENT 26
   * */
  for( j = n1 + nb; j < (n1 + 2*nb); j++ ){
    for( i = n1; i < nt; i++ ){
      ij = i + j*nt;
      ijh = i - n1 + (j + nc - n1)*n2;
      z[ij] = ftoc(h2[ijh],0.);
      for( l = 0; l < nb; l++ ){
	il = i - n1 + (l + nb + nc)*n2;
	lj = l + na + (j + na - n1 - nb)*nflt;
	z[ij] = csub(z[ij],fcm(g2[il],zn[lj]));
      }
    }
  }

  /****ELEMENT 17
   * */
  for( j = n1 + 2*nb; j < (n1 + 2*nb + nc); j++ ){
    for( i = 0; i < n1; i++ ){
      ijz = i + j*nt;
      ijg = i + (j + na - n1 - 2*nb)*n1;
      z[ijz] = ftoc(-g1[ijg],0.);
    }
  }

  /****ELEMENT 27
   * */
  for( j = n1 + 2*nb; j < (n1 + 2*nb + nc); j++ ){
    for( i = n1; i < nt; i++ ){
      ijz = i + j*nt;
      ijg = i - n1 + (j + nc - n1)*n2;
      z[ijz] = ftoc(g2[ijg],0.);
    }
  }

  /****ELEMENT 18
   * */
  for( j = nt - nc; j < nt; j++ ){
    for( i = 0; i < n1; i++ ){
      ijz = i + j*nt;
      ijg = i + (j + n1 - nt)*n1;
      z[ijz] = cfm(ex,g1[ijg]);
    }
  }

  /****ELEMENT 28
   * */
  for( j = nt - nc; j < nt; j++ ){
    for( i = n1; i < nt; i++ ){
      ijz = i + j*nt;
      ijg = i - n1 + (j + nb - nt + nc)*n2;
      z[ijz] = ftoc(-g2[ijg],0.);
    }
  }
  /*Finished Loading z--------------------------------*/
  if(quick == 1){
    zgetrf_(&nt,&nt,z,&nt,IPIV,&INFO);

    //  printf(" INFO = %d \n",INFO);
    //  printf("IPIV \n");

    perms = 0;
    for( i =0; i < nt; i++)
      if(i+1 != IPIV[i])
	perms++;

    perms %=2;
    if(perms)
      det = ftoc(-1.0,0.0);
    else
      det = ftoc(1.0,0.0);
    
    for(i = 0; i < nt; i++){
      ii = i*(nt+1);
      det = cmul(det,z[ii]);
    }
 
    //  printf("Z\n");
    //  cout_(nt,nt,z);

    //    zgecon_(NORM,&nt,z,&nt,&ANORM,&RCOND,WORK,RWORK,&INFO);

    printf("(%e,%e) (%e,%e) INFO=%d \n", g.re, g.im, det.re, det.im, INFO);
    // printf("RCON = %e\n", RCOND);
    //  exit(0);
  return(det);
  }
  else if( quick == 0 ){
    JOBVL[0] = 'N';
    JOBVR[0] = 'N';
    zgeev_(JOBVL, JOBVR, &nt, z, &nt, zeign, VL, &nt,
	   VR, &nt, WORK, &LWORK, RWORK, &INFO);
    IMIN = 0;
    RMIN = HUGE;
    ZMIN = ftoc(HUGE,HUGE);
    for(i = 0; i < nt; i++){
      rval = cabs1(zeign[i]);
      if(rval < RMIN){
	RMIN = rval;
	IMIN = i;
	ZMIN = zeign[i];
      }
    }
    printf("(%e,%e) (%e,%e) INFO=%d\n", g.re, g.im, ZMIN.re, ZMIN.im, INFO);
    return( ZMIN );
  }
  /* Now calculate the eigenvectors */
  else{
    /* Identify space harmonics leaking in superstrate and substrate */
    l = NLayers - 1;
    for(n = -nfl; n < nfl; n++){
      n_ = nfl+n;
      nl = n_;
      TopLeak[n_] = 0;
      if(h[nl].re < 0.0)
	TopLeak[n_] = 1;
      nl = n_ + l*nflt;
      BotLeak[n_] = 0;
      if(h[nl].re < 0.0 )
	BotLeak[n_] = 1;
    }
    /* Calculate the Matrix 1 norm. i.e. input NORM = '1' 
     * Used to calculate the condition number.  */
    //    ANORM = zlange_(NORM, &nt, &nt, z, &nt, RWORK);
    //  printf("ANORM = %e \n", ANORM);

    JOBVL[0] = 'N';
    JOBVR[0] = 'V';
    zgeev_(JOBVL, JOBVR, &nt, z, &nt, zeign, VL, &nt,
	   VR, &nt, WORK, &LWORK, RWORK, &INFO);
    printf("Eigenvectors Found.\n");
    OptWork = (int )WORK[0].re;
    printf("INFO = %d\tOptimal LWORK = %d\n", INFO, OptWork );
    if(LWORK < OptWork ){
      printf("LWORK = %d is less than the optical value.\n", LWORK );
      printf("Consider increasing its value.\n");
    }
    IMIN = 0;
    RMIN = HUGE;
    ZMIN = ftoc(HUGE,HUGE);
    for(i = 0; i < nt; i++){
      rval = cabs1(zeign[i]);
      if(rval < RMIN){
	RMIN = rval;
	IMIN = i;
	ZMIN = zeign[i];
      }
    }
    /* Transfer the zero eigenvector from the eigenvector matrix */
    for(i = 0; i < nt; i++)
      Evec[i] = VR[i+nt*IMIN];
    printf("Minimum Eigenvalue = (%e, %e)\n", ZMIN.re, ZMIN.im);
    /* The eigenvector has been normalized, i.e. Norm  = 1.
     * Now, normalize relative to value at the center of the
     * tooth where the field shoud be maximum.  Nevertheless,
     * this is an arbitrary condition. */
    icenter = (na-1)/2;  /* This requires na to be odd. */
    vcenter = Evec[icenter];
    for(i = 0; i < nt; i++)
      Evec[i] = cdiv(Evec[i],vcenter);
    //    cout_(nt,1,Evec);
    /* 
     * Remove the propagation terms from all vectors.
     * Vbot = (U1,V1)
     * Vtop = (U3,V3)
     * 
     * ---------------------------------------------
     * The vector locations:
     * N0    N1   N2   N3   N4   N5   N6   N7   N8
     * |-na--|-nc-|-na-|-nb-|-nc-|-nb-|-nc-|-nc-| 
     * | u1  | u2 | u3 | v1 | v2 | v3 | p2 | q2 |
     * u and v are fields in left and right boxes.
     * p and q are fluxs in left and right boxes.
     *-------------------------------------------------*/
    N0 = 0;
    N1 = N0+na;
    N2 = N1+nc;
    N3 = N2+na;
    N4 = N3+nb;
    N5 = N4+nc;
    N6 = N5+nb;
    N7 = N6+nc;
    N8 = N7+nc;
    /* Lift propagation. Ref to center of tooth */
    for(l = 0; l < nflt; l++)
      zgam[l] = zexp(cfm(g,(2*l+1-na)*dd2));
    /* u1 and u3 */ 
    for(l = 0; l < na; l++){
      Evec[l] = cmul(Evec[l],zgam[l]);
      Evec[l+N2] = cmul(Evec[l+N2],zgam[l]);
    }
    /* v1 and v3 */
    for(l = 0; l < nb; l++){
      Evec[l+N3] = cmul(Evec[l+N3],zgam[l+na]);
      Evec[l+N5] = cmul(Evec[l+N5],zgam[l+na]);
    }
    /* u2 and p2 */
    ex = zexp(cfm(g, na*dd2 ));
    for(l = 0; l < nc; l++){
      Evec[l+N1] = cmul(Evec[l+N1],ex);
      Evec[l+N6] = cmul(Evec[l+N6],ex);
    }
    /* v2 and q2 */
    ex = zexp(cfm(g,(na+2*nb)*dd2));
    for(l = 0; l < nc; l++){
      Evec[l+N4] = cmul(Evec[l+N4],ex);
      Evec[l+N7] = cmul(Evec[l+N7],ex);
    } /* End Lift of propagation. */

    /* Fill the top and bottom vectors */
    for(l = 0; l < na; l++){
      Vbot[l] = Evec[l];
      Vtop[l] = Evec[l+N2];
    }
    for(l = 0; l < nb; l++){
      Vbot[l+na] = Evec[l+N3];
      Vtop[l+na] = Evec[l+N5];
    } /* End Fill top and bottom vectors */

    /* Output the top and bottom vectors */
    printf(" \nThe Eigenvectors...(x + jy)\n");
    printf("   Vbot                              Vtop\n" );
    for(i = 0; i < nflt ; i++ ){
      printf(" %e  %e \t %e  %e \n",
	     Vbot[i].re, Vbot[i].im, Vtop[i].re, Vtop[i].im);
    }
    printf(" \nThe eigenvectors...(Mag and Phase)\n");
    printf("   Vbot                              Vtop\n" );
    for(i = 0; i < nflt ; i++) {
      Vbotmag = zabs(Vbot[i]);
      Vtopmag = zabs(Vtop[i]);
      Vbotpha = atan2(Vbot[i].im,Vbot[i].re);
      Vtoppha = atan2(Vtop[i].im,Vtop[i].re);
      printf("(%e, %e)\t(%e, %e)\n", 
	     Vbotmag, Vbotpha, Vtopmag, Vtoppha);
    }

    /* Compute Vtopp and Vbotp */
    convolve(nflt, Vtopp, tcn, Vtop);
    convolve(nflt, Vbotp, tcm, Vbot);
    
    /* Begin loading q1, p1, q2, p2 
     * Follow the Evec arrangement of vectors.
     * Propagation terms for fluxes must be calculated
     * from the Vtopp and Vbotp vectors.
     * -----------------------*/

    /*----1-----------------*/
    for(l = 0; l < na; l++){
      q1[l] = Evec[l];
      p1[l] = Vbotp[l];
    }
    /*----2-----------------*/
    for(l = 0; l < nc; l++){
      q1[l+na] = Evec[l+N1];
      p1[l+na] = Evec[l+N6];
    }
    /*----3-----------------*/
    for(l = 0; l < na; l++){
      q1[l+na+nc] = Evec[l+N2];
      p1[l+na+nc] = Vtopp[l];
    }/* Completed first three */
    /*----1-----------------*/
    for(l = 0; l < nb; l++){
      q2[l] = Evec[l+N3];
      p2[l] = Vbotp[l+na];
    }
    /*----2-----------------*/
    for(l = 0; l < nc; l++){
      q2[l+nb] = Evec[l+N4];
      p2[l+nb] = Evec[l+N7];
    }
    /*----3-----------------*/
    for(l = 0; l < nb; l++){
      q2[l+nb+nc] = Evec[l+N5];
      p2[l+nb+nc] = Vtopp[l+na];
    }/* Completed next three */

    for(l = 0; l < nc; l++){
      q2[l+2*nb+nc] =      q1[l+na];
      p2[l+2*nb+nc] = cneg(p1[l+na]);
    }
    for(l = 0; l < nc; l++){
      q1[l+2*na+nc] =      q2[l+nb];
      p1[l+2*na+nc] = cneg(p2[l+nb]);
    } /* Finished loading q1,p1,q2,p2 */
    /*
      cout_(n1,1,q1);
      cout_(n1,1,p1);
      cout_(n2,1,q2);
      cout_(n2,1,p2);
    */
    return( ZMIN );
  }
}
/*----------------------------------------------------------------------- */
void floquet( int Nt, int Nf, int Ntop, int Nbot, int NL,
	      complx *W, complx *WH, complx *Vt, complx *Vb,
	      complx *h, complx *c, complx *s, complx *t, 
	      complx *A, complx *B )
     /*-------------------------------------------- 
      * Nt Total floquet components
      * Nf- Floquet components left and right
      * Ntop - Top grating layer
      * Nbot - Bottom grating layer
      * W_n - exp( -j 2*Pi*n/N )
      * WH_n - exp( -j Pi*n/N )
      * Vt_n - Field at top of grating
      * Vb_n - Field at bottom of grating
      * h_n - Transverse wave numbers
      * A_n - coefficient of cosine
      * B_n - Coefficient of sine
      *-------------------------------------------*/
{
  int l, n, n_, nl, nn, nlm1;

  static complx J = {0.0,1.0};
  static double eps = 1.e-8;

  complx w;
  complx *v;
  complx Anl, Bnl;
  complx Ares, Bres;
  double dnflt;

  /* Now reference to the center of the left tooth */
  //icenter = (na-1)/2;
  //vcenter = cabs1(Vb[icenter]);
  //eps = 1.e-8;
  //eps = eps/vcenter;

  v = (complx *)malloc(Nt*sizeof(complx ));

  dnflt = (double )nflt;
  idft( Nt, Nf, W, Vb, v); /* Inverse Tr of bottom vectore */
  /* Calculate the psi_n (with adjusted phase referenced to
   * the center of the tooth). */
  l = Nbot;
  for( n = -Nf; n < Nf ; n++){
    n_ = Nf + n;
    nl = n_ + l*Nt;
    /* start phase adjuster, referenced to center of tooth */
    nn = (abs(n)*(na-1))%(2*nflt); 
    w = WH[nn];
    if( n < 0 )
      w = zconj(w);                
    /* end phase adjuster   */
    A[nl] = cmul(v[n_],w);
    B[nl] = cfm(cmul(cm[n_],A[nl]),-dnflt); /* normal is negative */
  }
  /* Space harmonics of bottom layers. 
   * Ares and Bres are residue values that are used
   * when the cosine terms c[nl] become vary large for high-order space
   * harmonics.   -------------------------------------------------- */
  for(l = Nbot+1; l < NL-1; l++) {
    for(n = -Nf; n < Nf ; n++) {
      n_ = Nf + n;
      nl = n_ + l*Nt;
      nlm1 = n_ + (l-1)*Nt;
      Anl = csub(A[nlm1],tdiv(J,B[nlm1],h[nl]));
      Ares = cneg(tdiv(B[nlm1],r[nl],h[nl]));
      Bnl = tmul(J,h[nl],Anl);
      Bres = tmul(A[nlm1],h[nl],r[nl]);
      //printf("(n,l), (Anl,Ares) ( %d, %d) (%e, %e)\n",
      //     n,l,cabs1(Anl),cabs1(Ares)); 
      /* If Anl is too small, its accuracy is lost */
      if(cabs1(Anl) < eps) {
	A[nl] = cmul(c[nl],Ares);
	B[nl] = cmul(c[nl],Bres);
      }
      else{
	A[nl] = cmul(c[nl],cadd(Anl,Ares));
	B[nl] = cmul(c[nl],cadd(Bnl,Bres));
      }
    }
  }
  /* Fill in the bottom layer */
  l = NL-1;
  for(n = -Nf; n < Nf ; n++) {
    n_ = Nf + n;
    nl = n_ + l*Nt;
    nlm1 = n_ + (l-1)*Nt;
    A[nl] = A[nlm1];
  }
  /*-------------------------------------------------------------
   * Do the same for the top layers 
   *-------------------------------------------------------------*/
  idft( Nt, Nf, W, Vt, v);
  l = Ntop-1;
  for( n = -Nf; n < Nf ; n++) {
    n_ = Nf + n;
    nl = n_ + l*Nt;
    /* start phase adjuster on top side of tooth */
    nn = (abs(n)*(na-1))%(2*nflt); 
    w = WH[nn];
    if(n < 0)
      w = zconj(w);
    /* end phase adjuster   */
    A[nl] = cmul(v[n_],w);
    B[nl] = cfm(cmul(cn[n_],A[nl]),dnflt); /* normal is positive */
  }
  /* Fill the remaining values of A and B by transformations*/

  for(l = Ntop-1; l < 0; l--) {
    for(n = -Nf; n < Nf ; n++) {
      n_ = Nf + n;
      nl = n_ + l*Nt;
      nlm1 = n_ + (l-1)*Nt;
      A[nlm1] = cadd(cmul(c[nl],A[nl]),tdiv(s[nl],B[nl],h[nl]));
      B[nlm1] = csub(cmul(c[nl],B[nl]),tmul(h[nl],s[nl],A[nl]));
    }
  }

  printf("A coefficients\n");
  cout_(Nt, NL-1, A);
  printf("B coefficients\n");
  cout_(Nt, NL-1, B);
  free(v);
  return;
}
/*---------------------------------------------------------------------*/
void calgs(int store,  int n,  int *ihaf, 
	   int *idub, double pecmax, complx *z1, complx *z2, complx *z3, 
	   double *c, double *dc)
     /*---------------------------------------------------------------
      * store - (=1) stores the root
      * n - the number of iterations
      * ihaf - number of consecutive halves
      * idub - number of consecutive doubles
      * pecmax - maximum percentage change
      * z1, z2, z3, guessed values for next root
      * c - iteration variable
      * dc - iteration variable step
      *--------------------------------------------------------------*/
{
  int Stored;
  int i, indx[10], nnow;
  double a, b, dal, dbe, dcmax, ddc, delta, eps, omd, opd, 
    peca, pecb, ta[10], tb[10], tc[10];

  eps = 1.e-4;
  dcmax = 0.004;
  Stored = n;

  if( n > 10 )
    nnow = 10;
  else
    nnow = n;
  
  a = (*z3).re;
  b = (*z3).im;
  if( nnow > 3 ) {
    peca = fabs( a - al[nnow-2] )/a;
    pecb = fabs( b - be[nnow-2] )/b;
    printf( "AL(NNOW-1) = %g \n", al[nnow-2] );
    printf( " A = %g \n", a );
    printf( "PECA = %g  PECMAX %g \n", peca, pecmax );
    printf( "PECB = %g  PECMAX %g \n", pecb, pecmax );
    if( (a > eps && al[nnow-2] > eps && peca > pecmax) 
	|| (pecb > 0.01) ) {
      printf( "PEC TOO BIG, INTERVAL HALVED...\n" );
      *ihaf += 1;
      *idub = 0;
      store = 0;
      *c -= *dc;
      *dc *= 0.5e0;
    }
  }

  /* FIRST STORE ROOT IN ARRAYS WHEN IDUB.NE.0. * */
  if( store ) {
    if( Stored > 10 ){
      for( i = 0; i < nnow-1; i++ ){
	cl[i] = cl[i + 1];
	al[i] = al[i + 1];
	be[i] = be[i + 1];
      }
    }
    cl[nnow-1] = *c;
    al[nnow-1] = a;
    be[nnow-1] = b;
    
    /* ORDER THE ROOTS * */
    if( nnow > 2 ){
      indexx( nnow, cl, indx );
      for( i = 0; i < nnow; i++ ){
	tc[i] = cl[i];
	ta[i] = al[i];
	tb[i] = be[i];
      }
      if( *dc > 0.0 ){
	for( i = 0; i < nnow; i++ ){
	  cl[i] = tc[indx[i]];
	  al[i] = ta[indx[i]];
	  be[i] = tb[indx[i]];
	}
      }
      else{
	for( i = 0; i < nnow; i++ ){
	  cl[nnow - 1 - i] = tc[indx[i]];
	  al[nnow - 1 - i] = ta[indx[i]];
	  be[nnow - 1 - i] = tb[indx[i]];
	}
      }

      *ihaf = 0;
      *idub += 1;
      if( *idub > 4 ){
	*dc *= 2.e0;
	*idub = 0;
	if( fabs( *dc ) >= dcmax ){
	  *dc = dcmax;
	  *idub = 0;
	   printf("MAX INTERVAL ATTAINED%g \n", *dc );
	}
	else{
	  printf( "INTERVAL DOUBLED%g \n", *dc );
	}
      }
    }
  }

  /*     ENTER HERE FOR NONSTORE OPERATION
   * */
  *c += *dc;
  if( nnow < 3 ){
    delta = 1.e-5;
    omd = 1.e0 - delta;
    opd = 1.e0 + delta;
    *z1 = cadd(cmul(ftoc( omd, delta ),*z3),ftoc(delta,0.));
    *z2 = csub(cmul(ftoc( opd, delta ),*z3),ftoc(delta,0.));
  }
  else{
    *z1 = ftoc( al[nnow-1], be[nnow-1] );
    dal = al[nnow-1] - al[nnow - 2];
    dbe = be[nnow-1] - be[nnow - 2];
    ddc = cl[nnow-1] - cl[nnow - 2];
    a = al[nnow-1] + dal**dc/ddc;
    b = be[nnow-1] + dbe**dc/ddc;
    *z3 = ftoc( a, b );
    a = al[nnow-1] + 2.e0*dal**dc/ddc;
    b = be[nnow-1] + 2.e0*dbe**dc/ddc;
    *z2 = ftoc( a, b );
  }
  return;
} /* end of function */







