/* Copyright Jerome K. Butler             1992-2002
 *****************************************************/
typedef struct{
  double re;
  double im;
} complx;


static const complx ZZERO = {0.0, 0.0};
static const complx ZONE = {1.0, 0.0};
static const complx ZJ = {0.0, 1.0};

/* A complx 2x2 matrix */
typedef struct{
  complx A;
  complx B;
  complx C;
  complx D;
} Zmatrix;

typedef struct{
  complx a;
  complx b;
} Zvector;

/* Structure of Newton search */
typedef struct{
  double R; /* Reciprocal Matrix Condition */
  complx z;
} MatNew;

typedef struct{   /* Geometry ****************************/
  int L;          /* number of layers */
  int GL;         /* reference layer (bottom of grating */
  int GLF;        /* Layer with the fill material */
  int NLG;        /* Number of layers in grating region */
  double k;       /* free-space wavenumber */
  int nX;         /* number of x plot points*/
  /* arrays */
  char *C;        /* Layer Comments */
  int *LN;        /* Layer number identifier */
  int *LF;        /* Boolean: Layer Flat */
  int *Cells;     /* Number of repeated unit cells */
  int *CellLayers;/* Number of layers in each unit cell */
  double *w;      /* Layer thicknesses in units of lambda */
  double *d;      /* Normalized layer thicknesses */
  double *nl;     /* Left refractive index */
  double *nr;     /* Right refractive index */
  double *al;     /* wave absorption /cm */
  double *ext;    /* extinction coefficient */
  complx *Kl;     /* initial dielectric constants of layer */
  complx *Kr;     /* final dielectric constants of layer */
  double *bcept;  /* Intercept point on kappa axis */
  double *slp;    /* slope of kappa line */
  double *am;     /*                     */
  double *tm;     /*                     */
  double *xpts;   /* layer boundaries */
  double *xrpts;  /* reverse layer boundaries */
  double *Xf;     /* field plot points Xf[nX] */
} Geometry;

typedef struct{    /* Unloaded *****************************/
  int L;           /* number of layers */
  int GL;          /* Grating layer reference only */
  complx g;        /* Propagation Constant \gamma */
  /* arrays */
  complx *h;
  complx *c;
  complx *s;
  complx *t;
  complx *q;       /* Field values */
  complx *p;       /* Flux values */
  double *Gamma;       /* unloaded confinement */
} Unloaded;

typedef struct{/* Green *************************************/
  int Boxes;
  int *Boxpts;
  int *BoxXdim;
  int *BoxZdim;
  int na;
  int nb;
  int nc;
  int n1;  /* 2*(na + nc ) */
  int n2;  /* 2*(nb + nc ) */
  int nt;  /* n1 + n2 SAME AS EIGEN.N */
  int nts; /* nt*nt SAME AS EIGEN.N2 */
  /* Updated parameters with each itt */
  double CL;  /* Normalized grating period */
  double K;   /* Normalized grating wavenumber */
  double del;
  double delc;
  double dd2;
  double tauz;
  /* arrays */
  int *Ix; /* Panel centers along x Ix[i] * delc/2 */
  int *Iz; /* Panel centers along z Iz[i] * del/2  */
  int *Ir; /* Ir[ij] = Repeat panel */
  int *ix1;
  int *ix2;
  int *iz1;
  int *iz2;
  int *ir1;
  int *ir2;
  double *G;
  double *H;
  double *g1;
  double *g2;
  double *h1;
  double *h2;
} Green;

typedef struct{ /* DFT ***************************************/
  int nfl;
  int NF;
  complx *W;
  complx *PF;
  complx *TN;
  complx *WH;
  complx *phin;
  complx *Phil;
  complx *FpPhil;
} DFT;

typedef struct{ /* Loaded **************************************/
  int L;      /* Number of layers */
  int GL;     /* Grating layer */
  int nfl;    /* NF > 1 */
  int NF;     /* Number of Space harmonics */
  int LBG;    /* number of layers below the grating */
  int nrts;   /* number of Fncn fuction calls */
  complx g;   /* Same as gn(n=0) */
  /*----------------------------------------------------------*/
  /* Below are dimensioned by NF */
  char *TopLeak;
  char *BotLeak;
  complx *gn;
  complx *zgam;
  complx *cm;
  complx *cmp;
  complx *tcm;
  complx *tcmp;
  complx *cn;
  complx *cnp;
  complx *tcn;
  complx *tcnp;
  /* Below are dimensioned by NF*NF */
  complx *zm;
  complx *zmp;
  complx *zn;
  complx *znp;
  complx *U;
  complx *UP;
  complx *O;
  complx *Emgz;
  /* Below are dimensioned by L*NF (Floquet n changes fastest)*/
  char *illCond;
  complx *h;
  complx *hp;
  complx *hd;
  complx *e;
  complx *c;
  complx *cp;
  complx *s;
  complx *sp;
  complx *soh;
  complx *sth;
  complx *sohp;
  complx *sthp;
  complx *t;
  complx *toh;
  complx *tth;
  complx *tohp;
  complx *tthp;
  complx *q;
  complx *p;
  complx *BK;   /* bra-ket elements */
  double *Gamma;/*Layer confinemnt */
  Zvector *Tp;  /* dB/dg and dC/dg elements of T' */
  Zmatrix *TT;  /* Transfer matrix of top layers */
  Zmatrix *TB;  /* Transfer matrix of bottom layers */
  Zmatrix *TTP; /* Derivative of top transfer matrix */
  Zmatrix *TBP; /* Derivative of bottom transfer matrix */
  Zmatrix *BR;  /* Bottom right partial products */
  Zmatrix *BL;  /* Bottom left partial products */
  Zmatrix *TR;  /* Top right recurrence matrices */
  Zmatrix *TL;  /* Top left recurrence matrices */
} Loaded;

typedef struct{ /* Eigen **************************************/
  int NU;  /* Unloaded matrix dimension 2*(GL-1) */
  int NL; /* OLD nt Loaded dimension */
  int NLS; /* Super matrix dimension */
  int IG; /* Beginning row for grating matrix */
  int JG; /* Beginning col for grating matirx */
  int IB; /* Bottom of the grating */
  int JB; /* Top of the grating */
  int quick; /* Boolean variable for solving eigenvalues */
  int LWORK;
  int INFO;
  double RCOND;
  /* arrays */
  complx *A;  /* Master matrix: L x L,  the old z */
  complx *B; /* Derivative w.r.t \gamma the old zp */
  complx *VR; /* Right eigenvectors */
  complx *VL; /* Lest eigenvectors */
  complx *AF;
  complx *X;
  complx *WORK;
  complx *W;  /* Eigenvalues on diagonal (ikd zeign) */
  complx *Evec;    /* Boundary Element Eigenvector as below*/
  double *RWORK;
  double *BERR;
  double *C;
  double *FERR;
  double *R;
  int    *IPIV;
  char   *JOBVL;
  char   *JOBVR;
  char   *EQUED;
  char   *FACT;
  char   *TRANS;
} Eigen;

typedef struct{ /* GratFields *********************************/
  int na;
  int nb;
  int nc;
  complx *Vtop;    /* Ref to na/2; Propagation removed in fncn */
  complx *Vtopp;
  complx *Vbot;    /* Ref to na/2; Propagation removed in fncn */
  complx *Vbotp;
  /* Counting reversed on top and left sides of boxes */
  complx *q1;      /* Field points around 1st box */
  complx *p1;      /* Out flux around 1st box */
  complx *q2;      /* Field points around 2nd box */
  complx *p2;      /* Out flux around 2nd box */
  complx *Fgrat;   /* F nc*(na+nb) = \sum \psi_n(x) exp(-jnKz) */
  complx *FPgrat;  /* FP  = \sum \n psi_n(x) exp(-jnKz) */
  double *Maggrat; /* Man of Fgrat */
  double *Phagrat; /* phase of Fgrat */
  double *xgrat;   /* x locations in the gratins nc*(na+nb) */
  double *zgrat;   /* z locations in the grating nc*(na+nb) */
  /* Starting Spline variables */
  double *xspl;    /* x locations in the grating nc+2 */
  double *psre;    /* psi.re at the xspl points */
  double *psim;    /* psi.im at the xspl points */
  double *bre;     /* |---------------------------------------*/
  double *cre;     /* | Real Spline coefficients of grat fld. */
  double *dre;     /* |---------------------------------------*/
  double *bim;     /* |--------------------------------------------*/
  double *cim;     /* | Imaginary Spline coefficients of grat fld. */
  double *dim;     /* |--------------------------------------------*/
  double *psrep;   /* F'.re at the xspl points */
  double *psimp;   /* F'.im at the xspl points */
  double *brep;    /* |---------------------------------------*/
  double *crep;    /* | Real Spline coefficients of F'        */
  double *bimp;    /* |--------------------------------------------*/
  double *cimp;    /* | Imaginary Spline coefficients of F'        */
}GratFields;
  

/* complx function prototypes */

complx cneg( complx );   /* negation of a dp complx no. */
/*ADDITION*/
complx cadd( complx, complx ); /* sum of two complx nos. */
complx tadd( complx, complx, complx ); /* sum of two complx nos. */
complx cfadd( complx, double );/* add a complx to a double. */
complx fcadd( double, complx );/* add a double to a complx. */
/*SUBTRACTION*/
complx csub( complx, complx ); /* difference of two complx nos. */
complx tsub( complx, complx, complx); /* A - B - C */
complx cfsub( complx, double );/* sub a double from a complx. */
complx fcsub( double, complx );/* sub a complx from a double. */
/*MULTIPLICATION*/
complx csq( complx );          /*Square of a complex number */
complx zcub( complx );         /* Cube of complex number */
complx cmul( complx, complx ); /* product of two complx nos. */
complx fcm( double, complx );  /* product of double and  complx nos. */
complx cfm( complx, double );  /* product of complx and double */
complx tmul( complx, complx, complx );   /* product of three complx nos. */
complx fmul( complx, complx, complx, complx );/* product of four complx nos.*/
complx cdm( complx, complx, complx, double );/* product of four complx nos.*/
/*DIVISION*/
complx cdd(complx, double );   /* complx divided by double. */
complx dcd(double, complx );   /* double divided by complx. */
complx cdiv( complx, complx ); /* quotient of two complx nos. */
complx tdiv( complx, complx, complx );   /* a*b/c */
complx fdiv( complx, complx, complx, complx );   /* a*b*c/d */
complx zinv( complx ); /* inverse of a complx number. */

complx dconjg( complx ); /* complx conjugate */
/* Muller's method - zeros of a complx function */
complx croot(complx (* )(complx ),complx *, complx *, int *, int , int *);
/* Newton's method - zeros of a complx function */
complx znewt(MatNew (* )(complx ), complx, double, int *, int , int *);
complx Znewton(MatNew(*fcn)(complx , Geometry , Unloaded , Eigen ), 
	       Geometry G, Unloaded U, Eigen E,
	       complx z, double eps, int *i, int m, int *root);


int ceq( complx, complx );   /* return TRUE if complx nos. equal */

   /* complx Library Functions */
complx ftoc(double,double); /* convert floats to dp complx */
complx csqrt( complx );  /* (positive) sq root of dp complx no. */
complx zconj( complx );  /* dp complx conjugate */
complx zexp( complx );   /* dp complx exponential */
complx zlog( complx );   /* dp complx natural log */
complx zsqrt( complx );  /* branch cut in third quadrant. */
complx ztan( complx );   /* tangent of dp complx no. */
complx zunounp1(int , complx ); /* Ratio of Cheb polys. */
complx casin( complx );  /* sine of dp complx no. */
complx cacos( complx );  /* sine of dp complx no. */
complx csin( complx );   /* sine of dp complx no. */
complx ccos( complx );   /* cosine of dp complx no. */
complx btan( complx );   /* tangent residue with a big arg. */
complx czeta( complx );  /* czeta(z) = tan(z) - j  */

double zabs( complx );   /* dp complx absolute value */
double zabs2( complx );   /* dp complx absolute value squared.*/
double cabs1( complx );  /* dp sum of |real| and |imag| parts. */

void cguesses(complx *);
void zcheb( int , complx , complx *, complx * ); /*Un(z) chebyshev */

Zvector zvm(Zmatrix , Zvector); /* Multiply complx matrix and complx vec*/
Zmatrix zadd(Zmatrix, Zmatrix); /* Add 2 matrices */
Zmatrix zmm(Zmatrix , Zmatrix); /* Multiply two complex matrices */
Zmatrix ztm(Zmatrix , Zvector , Zmatrix); /* zvector is cross diag matrix */

/* Allocate complx array check */
void ComplexAllocated(complx *C, char *Cc);
/* Allocate a complx array and also check allocation */
complx *AllocateComplx(int N, char *Cc);
Zvector *AllocateZvector(int N, char *Cc);
Zmatrix *AllocateZmatrix(int N, char *Cc);
/* Complex Integrations along the real axis */

complx zromb(complx (*eval)(double ), double a, double b, double eps);

