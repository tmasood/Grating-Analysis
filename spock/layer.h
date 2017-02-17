#ifndef LAYER_INCLUDE
#define LAYER_INCLUDE

#include "local_complex.h"

struct LAYER
{
  int matsys;
  float xperc;
  float yperc;
  double nreal;
  double nloss;
  double per;
  double pei;
  double pmr;
  double pmi;
  double tl; /* layer thickness */
  double tr; /* radian thickness */
  double bl; /* length at boundry */
  dcomplex wx; /* transverse propagation coeff for each layer */
  dcomplex yn; /* normalized wave admittance(TE) / impedance(TM) factor */
  dcomplex yx; 
  dcomplex yz;
  dcomplex kappa;
  dcomplex ph;
  dcomplex ph0; /* complex phase thickness */
  dcomplex qn; /* normalized propagation constant k**2 or n**2 */
  dcomplex rn; /* sqrt(qn) = n */
  dcomplex cl0; /* matrix entries for each layer */
  dcomplex cl1; /* matrix entries for each layer */
  dcomplex cl2; /* matrix entries for each layer */
  dcomplex fy; /* boundry field value */
  dcomplex fz; /* boundry field value */
  dcomplex fym[MAXLOOP]; /* store boundry field value */
  dcomplex fzm[MAXLOOP]; /* store boundry field value */
  dcomplex gy; /* boundry field value */
  dcomplex gz; /* boundry field value */
  dcomplex px;
  double xl;
  double ra[MAXMODES];
  struct LAYER *nextptr;
  struct LAYER *prevptr;
};

struct GLAYER
{
  int gradedmatsys;
  float xperci;
  float yperci;
  float xpercf;
  float ypercf;
  double nreali;
  double nrealf;
  double nloss;
  double gtl;
  int nslc;
};

#endif
