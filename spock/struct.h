#include "local_complex.h"

#ifndef STRUCT_INCLUDE
#define STRUCT_INCLUDE

struct STRUCT
{
  float wvl;
  dcomplex kc; /* frequency factor */
  double k0;
  dcomplex kf; /* relative frequency wr */
  int mn; /* total number of modes to be searched */
  int mny; /* number of y modes */
  int ny;
  int mm; /* mode index  */
  int mo;
  int l1; /* each side outer boundry */
  int l2; /* outer boundry */
  double qzr; /* real part of neff**2(beta**2 used interchangeably) */
  double qzi; /* imaginary part of neff**2 */
  dcomplex qz; /* neff**2 */
  dcomplex wz; /* neff */
  dcomplex wzm[MAXMODES]; /* neff */
  dcomplex qx; /* transverse neff**2 */
  dcomplex phr; /* real part of phase integral */
  dcomplex phi; /* imaginary part of phase integral */
  dcomplex detc[2]; /* determinant and antideterminant */
  dcomplex dets;
  double atqz; /* polar phase angle of qz */
  int totallayers;
  double fwhpn;
  double fwhpf;
  int lzcnt[MAXSIMLOOP];
  int loopzv;
  double zfinv[MAXSIMLOOP][MAXSIMLOOP];
  double zinc[MAXSIMLOOP][MAXSIMLOOP];
  int ilz[MAXSIMLOOP][MAXSIMLOOP];
};

#endif
