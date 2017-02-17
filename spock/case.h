#ifndef CASE_INCLUDE
#define CASE_INCLUDE

#include "util.h"
#include "local_complex.h"

struct CASE
{
  int kase;
  int ksub;
  int last;
  int ktuo;
  float eps1;
  float eps2;
  float gameps;
  float qzmr;
  float qzmi;
  float intqzmr; /* initial root guess (real part)*/
  float intqzmi; /* initial root guess (imaginary part)*/
  float pmfr;
  float pmdm;
  dcomplex qzm[MAXMODES]; /* initial guesses for all modes */
  dcomplex z0[MAXMODES]; /* input used for guesses */
  dcomplex phm[MAXMODES]; /* normalized phase integral */
  float fm[MAXMODES]; /* absolute magnitude of function at roots */
  int it[MAXMODES]; /* number of iterations completed for each root */
  int km[MAXMODES]; /* saved kr */
  int kr[MAXMODES]; /* convergence indicator (root status and quality) */
  int kmdo;
  float qznr; /* initial guess of real part of neff**2 */
  float qzni; /* initial guess of imaginary part of neff**2 */
  float printf;
  float initgs;
  float autoqw;
  float nfplt;
  float ffplt;
  float aqz;
  float dxin;
  double dtheta;
  double thetax;
  int peropt;
  int lxyopt;
  int kdofy;
  int koufy;
  int kgczy;
  int kgssy;
  int il; /* iteration limit */
  int kgss; /* type of initial guess */
  int kgcz; /* how guess is used in complex root search */
  int kdof; /* calculate field */
  int kouf; /* output field */
  int keif;
  int pfield;
  int ntheta;
};

#endif
