/************************************************************************
 *                                                                      *
 * File: cvkxb.c                                                        *
 * Programmers: Scott D. Cohen and Alan C. Hindmarsh @ LLNL             *
 * Version of 23 March 2000                                             *
 *----------------------------------------------------------------------*
 * Example problem.                                                     *
 * An ODE system is generated from the following 2-species diurnal      *
 * kinetics advection-diffusion PDE system in 2 space dimensions:       *
 *                                                                      *
 * dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)    *
 *                 + Ri(c1,c2,t)      for i = 1,2,   where              *
 *   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,       *
 *   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,                    *
 *   Kv(z) = Kv0*exp(z/5) ,                                             *
 * Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)        *
 * vary diurnally.   The problem is posed on the square                 *
 *   0 <= x <= 20,    30 <= z <= 50   (all in km),                      *
 * with homogeneous Neumann boundary conditions, and for time t in      *
 *   0 <= t <= 86400 sec (1 day).                                       *
 * The PDE system is treated by central differences on a uniform        *
 * 10 x 10 mesh, with simple polynomial initial profiles.               *
 * The problem is solved with CVODE, with the BDF/GMRES method (i.e.    *
 * using the CVSPGMR linear solver) and a banded preconditioner,        *
 * generated by difference quotients, using the module CVBANDPRE.       *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "llnltyps.h"  /* definitions of real, integer, TRUE, FALSE       */
#include "cvode.h"     /* main CVODE header file                          */
#include "iterativ.h"  /* contains the enum for types of preconditioning  */
#include "cvspgmr.h"   /* use CVSPGMR linear solver each internal step    */
#include "cvbandpre.h" /* band preconditioner function prototypes         */
#include "nvector.h"   /* definitions of type N_Vector, macro N_VDATA     */
#include "llnlmath.h"  /* contains SQR macro                              */


/* Problem Constants */

#define NUM_SPECIES  2             /* number of species         */
#define KH           4.0e-6        /* horizontal diffusivity Kh */
#define VEL          0.001         /* advection velocity V      */
#define KV0          1.0e-8        /* coefficient in Kv(z)      */
#define Q1           1.63e-16      /* coefficients q1, q2, c3   */ 
#define Q2           4.66e-16
#define C3           3.7e16
#define A3           22.62         /* coefficient in expression for q3(t) */
#define A4           7.601         /* coefficient in expression for q4(t) */
#define C1_SCALE     1.0e6         /* coefficients in initial profiles    */
#define C2_SCALE     1.0e12

#define T0           0.0           /* initial time */
#define NOUT         12            /* number of output times */
#define TWOHR        7200.0        /* number of seconds in two hours  */
#define HALFDAY      4.32e4        /* number of seconds in a half day */
#define PI       3.1415926535898   /* pi */ 

#define XMIN          0.0          /* grid boundaries in x  */
#define XMAX         20.0           
#define ZMIN         30.0          /* grid boundaries in z  */
#define ZMAX         50.0
#define XMID         10.0          /* grid midpoints in x,z */          
#define ZMID         40.0

#define MX           10             /* MX = number of x mesh points */
#define MZ           10             /* MZ = number of z mesh points */
#define NSMX         20             /* NSMX = NUM_SPECIES*MX */
#define MM           (MX*MZ)        /* MM = MX*MZ */

/* CVodeMalloc Constants */

#define RTOL    1.0e-5            /* scalar relative tolerance */
#define FLOOR   100.0             /* value of C1 or C2 at which tolerances */
                                  /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)      /* scalar absolute tolerance */
#define NEQ     (NUM_SPECIES*MM)  /* NEQ = number of equations */


/* User-defined vector and matrix accessor macros: IJKth, IJth */

/* IJKth is defined in order to isolate the translation from the
   mathematical 3-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. IJth is defined in order to
   write code which indexes into small dense matrices with a (row,column)
   pair, where 1 <= row, column <= NUM_SPECIES.   
   
   IJKth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
   0 <= j <= MX-1, 0 <= k <= MZ-1. The vdata array is obtained via
   the macro call vdata = N_VDATA(v), where v is an N_Vector. 
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata.

   IJth(a,i,j) references the (i,j)th entry of the small matrix real **a,
   where 1 <= i,j <= NUM_SPECIES. The small matrix routines in dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. */


#define IJKth(vdata,i,j,k) (vdata[i-1 + (j)*NUM_SPECIES + (k)*NSMX])
#define IJth(a,i,j)        (a[j-1][i-1])


/* Type : UserData 
   contains preconditioner blocks, pivot arrays, and problem constants */

typedef struct {
  real **P[MX][MZ], **Jbd[MX][MZ];
  integer *pivot[MX][MZ];
  real q4, om, dx, dz, hdco, haco, vdco;
} *UserData;


/* Private Helper Functions */

static void InitUserData(UserData data);
static void SetInitialProfiles(N_Vector y, real dx, real dz);
static void PrintOutput(long int iopt[], real ropt[], N_Vector y, real t);
static void PrintFinalStats(long int iopt[]);

/* Function Called by the CVODE Solver */

static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);


/***************************** Main Program ******************************/

main()
{
  real abstol, reltol, t, tout, ropt[OPT_SIZE];
  long int iopt[OPT_SIZE];
  N_Vector y;
  UserData data;
  CVBandPreData bpdata;
  void *cvode_mem;
  int ml, mu, iout, flag;

  /* Allocate memory, and set problem data, initial values, tolerances */ 

  y = N_VNew(NEQ, NULL);
  data = (UserData) malloc(sizeof *data);
  InitUserData(data);
  SetInitialProfiles(y, data->dx, data->dz);
  abstol=ATOL; reltol=RTOL;

  /* Call CVBandPreAlloc to initialize band preconditioner */
  ml = mu = 2;
  bpdata = CVBandPreAlloc (NEQ, f, data, mu, ml);

  /* Call CVodeMalloc to initialize CVODE: 
     NEQ     is the problem size = number of equations
     f       is the user's right hand side function in y'=f(t,y)
     T0      is the initial time
     y       is the initial dependent variable vector
     BDF     specifies the Backward Differentiation Formula
     NEWTON  specifies a Newton iteration
     SS      specifies scalar relative and absolute tolerances
     &reltol and &abstol are pointers to the scalar tolerances
     data    is the pointer to the user-defined block of coefficients
     FALSE   indicates there are no optional inputs in iopt and ropt
     iopt    and ropt arrays communicate optional integer and real input/output

     A pointer to CVODE problem memory is returned and stored in cvode_mem.  */

  cvode_mem = CVodeMalloc(NEQ, f, T0, y, BDF, NEWTON, SS, &reltol,
                          &abstol, data, NULL, FALSE, iopt, ropt, NULL);
  if (cvode_mem == NULL) { printf("CVodeMalloc failed."); return(1); }

  /* Call CVSpgmr to specify the CVODE linear solver CVSPGMR with
     left preconditioning, modified Gram-Schmidt orthogonalization,
     default values for the maximum Krylov dimension maxl and the tolerance
     parameter delt, preconditioner setup and solve routines CVBandPrecond
     and CVBandPSolve, and the pointer to the user-defined block data.     */

  CVSpgmr(cvode_mem, LEFT, MODIFIED_GS, 0, 0.0, CVBandPrecond, CVBandPSolve,
          bpdata);

  printf("2-species diurnal advection-diffusion problem, %d by %d mesh\n", MX, MZ);
  printf("SPGMR solver; band preconditioner on left; mu, ml = %d  %d\n\n", mu, ml);

  /* In loop over output points, call CVode, print results, test for error */

  for (iout = 1, tout = TWOHR; iout <= NOUT; iout++, tout += TWOHR) {
    flag = CVode(cvode_mem, tout, y, &t, NORMAL);
    PrintOutput(iopt, ropt, y, t);
    if (flag != SUCCESS) { printf("CVode failed, flag=%d.\n", flag); break; }
  }

  /* Free memory and print final statistics */  

  N_VFree(y);
  CVBandPreFree(bpdata);
  CVodeFree(cvode_mem);
  PrintFinalStats(iopt);
  return(0);
}


/*********************** Private Helper Functions ************************/

/* Load problem constants in data */

static void InitUserData(UserData data)
{
  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/(MX-1);
  data->dz = (ZMAX-ZMIN)/(MZ-1);
  data->hdco = KH/SQR(data->dx);
  data->haco = VEL/(2.0*data->dx);
  data->vdco = (1.0/SQR(data->dz))*KV0;
}

/* Set initial conditions in y */

static void SetInitialProfiles(N_Vector y, real dx, real dz)
{
  int jx, jz;
  real x, z, cx, cz;
  real *ydata;

  /* Set pointer to data array in vector y. */

  ydata = N_VDATA(y);

  /* Load initial profiles of c1 and c2 into y vector */

  for (jz = 0; jz < MZ; jz++) {
    z = ZMIN + jz*dz;
    cz = SQR(0.1*(z - ZMID));
    cz = 1.0 - cz + 0.5*SQR(cz);
    for (jx = 0; jx < MX; jx++) {
      x = XMIN + jx*dx;
      cx = SQR(0.1*(x - XMID));
      cx = 1.0 - cx + 0.5*SQR(cx);
      IJKth(ydata,1,jx,jz) = C1_SCALE*cx*cz; 
      IJKth(ydata,2,jx,jz) = C2_SCALE*cx*cz;
    }
  }
}

/* Print current t, step count, order, stepsize, and sampled c1,c2 values */

static void PrintOutput(long int iopt[], real ropt[], N_Vector y, real t)
{
  real *ydata;

  ydata = N_VDATA(y);

  printf("t = %.2e   no. steps = %d   order = %d   stepsize = %.2e\n",
         t, iopt[NST], iopt[QU], ropt[HU]);
  printf("c1 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n",
         IJKth(ydata,1,0,0), IJKth(ydata,1,4,4), IJKth(ydata,1,9,9));
  printf("c2 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n\n",
         IJKth(ydata,2,0,0), IJKth(ydata,2,4,4), IJKth(ydata,2,9,9));
}

/* Print final statistics contained in iopt */

static void PrintFinalStats(long int iopt[])
{
  printf("\nFinal Statistics.. \n\n");
  printf("lenrw   = %5ld    leniw = %5ld\n", iopt[LENRW], iopt[LENIW]);
  printf("llrw    = %5ld    lliw  = %5ld\n", iopt[SPGMR_LRW], iopt[SPGMR_LIW]);
  printf("nst     = %5ld    nfe   = %5ld\n", iopt[NST], iopt[NFE]);
  printf("nni     = %5ld    nli   = %5ld\n", iopt[NNI], iopt[SPGMR_NLI]);
  printf("nsetups = %5ld    netf  = %5ld\n", iopt[NSETUPS], iopt[NETF]);
  printf("npe     = %5ld    nps   = %5ld\n", iopt[SPGMR_NPE], iopt[SPGMR_NPS]);
  printf("ncfn    = %5ld    ncfl  = %5ld\n \n", iopt[NCFN], iopt[SPGMR_NCFL]);
}


/***************** Functions Called by the CVODE Solver ******************/

/* f routine. Compute f(t,y). */

static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data)
{
  real q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt;
  real c1rt, c2rt, czdn, czup, hord1, hord2, horad1, horad2;
  real qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, zdn, zup;
  real q4coef, delz, verdco, hordco, horaco;
  real *ydata, *dydata;
  int jx, jz, idn, iup, ileft, iright;
  UserData data;

  data = (UserData) f_data;
  ydata = N_VDATA(y);
  dydata = N_VDATA(ydot);

  /* Set diurnal rate coefficients. */

  s = sin(data->om*t);
  if (s > 0.0) {
    q3 = exp(-A3/s);
    data->q4 = exp(-A4/s);
  } else {
    q3 = 0.0;
    data->q4 = 0.0;
  }

  /* Make local copies of problem variables, for efficiency. */

  q4coef = data->q4;
  delz = data->dz;
  verdco = data->vdco;
  hordco  = data->hdco;
  horaco  = data->haco;

  /* Loop over all grid points. */

  for (jz = 0; jz < MZ; jz++) {

    /* Set vertical diffusion coefficients at jz +- 1/2 */

    zdn = ZMIN + (jz - .5)*delz;
    zup = zdn + delz;
    czdn = verdco*exp(0.2*zdn);
    czup = verdco*exp(0.2*zup);
    idn = (jz == 0) ? 1 : -1;
    iup = (jz == MZ-1) ? -1 : 1;
    for (jx = 0; jx < MX; jx++) {

      /* Extract c1 and c2, and set kinetic rate terms. */

      c1 = IJKth(ydata,1,jx,jz); 
      c2 = IJKth(ydata,2,jx,jz);
      qq1 = Q1*c1*C3;
      qq2 = Q2*c1*c2;
      qq3 = q3*C3;
      qq4 = q4coef*c2;
      rkin1 = -qq1 - qq2 + 2.0*qq3 + qq4;
      rkin2 = qq1 - qq2 - qq4;

      /* Set vertical diffusion terms. */

      c1dn = IJKth(ydata,1,jx,jz+idn);
      c2dn = IJKth(ydata,2,jx,jz+idn);
      c1up = IJKth(ydata,1,jx,jz+iup);
      c2up = IJKth(ydata,2,jx,jz+iup);
      vertd1 = czup*(c1up - c1) - czdn*(c1 - c1dn);
      vertd2 = czup*(c2up - c2) - czdn*(c2 - c2dn);

      /* Set horizontal diffusion and advection terms. */

      ileft = (jx == 0) ? 1 : -1;
      iright =(jx == MX-1) ? -1 : 1;
      c1lt = IJKth(ydata,1,jx+ileft,jz); 
      c2lt = IJKth(ydata,2,jx+ileft,jz);
      c1rt = IJKth(ydata,1,jx+iright,jz);
      c2rt = IJKth(ydata,2,jx+iright,jz);
      hord1 = hordco*(c1rt - 2.0*c1 + c1lt);
      hord2 = hordco*(c2rt - 2.0*c2 + c2lt);
      horad1 = horaco*(c1rt - c1lt);
      horad2 = horaco*(c2rt - c2lt);

      /* Load all terms into ydot. */

      IJKth(dydata, 1, jx, jz) = vertd1 + hord1 + horad1 + rkin1; 
      IJKth(dydata, 2, jx, jz) = vertd2 + hord2 + horad2 + rkin2;
    }
  }
}
