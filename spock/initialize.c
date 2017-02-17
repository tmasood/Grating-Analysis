#include <stdio.h>
#include "util.h"
#include "case.h"
#include "struct.h"
#include "modcon.h"

int initialize(struct CASE *caseptr, struct STRUCT *structptr,
	   struct MODCON *modconptr)
{
  int m, i;
  dcomplex CZERO, CONE;

  __real__ CZERO = 0.0;
  __imag__ CZERO = 0.0;
  __real__ CONE = 1.0;
  __imag__ CONE = 0.0;

  structptr->mn = 1;
  structptr->mny = 1;
  structptr->ny = 0;

  structptr->phr = CZERO;
  structptr->phi = CZERO;

  for (m = 0; m < structptr->mn; m++)
    {
      caseptr->km[m] = 5;
      __real__ caseptr->qzm[m] = 10.65;
      __imag__ caseptr->qzm[m] = 0.0;
    }
  structptr->mo = 0;
  caseptr->qznr = 0.9;
  caseptr->qzni = -1.0e-4;
  caseptr->kgss = 1;
  caseptr->kgcz = 6;
  caseptr->eps1 = 1.0e-8;
  caseptr->eps2 = 1.0e-8;
  caseptr->il = 30;
  caseptr->lxyopt = 0;
  caseptr->peropt = 1;
  caseptr->kdofy = 1;
  caseptr->koufy = 3;
  caseptr->kgczy = 3;
  caseptr->kgssy = 1;
  caseptr->kase = 1;
  caseptr->ksub = 0;
  caseptr->last = 0;
  caseptr->ktuo = 4;
  caseptr->keif = 1;
  caseptr->pmfr = 0.3;
  caseptr->pmdm = 1.0;
  caseptr->kmdo = 5;
  caseptr->dtheta = 1.0;
  caseptr->thetax = 90.0;
  caseptr->dxin = 0.0;

  __real__ structptr->kc = 1.0;
  __imag__ structptr->kc = 0.0;
  __real__ structptr->kf = 1.0;
  __imag__ structptr->kf = 0.0;

  modconptr->kpol = 0;
  modconptr->kbc1 = 1;
  modconptr->kbc2 = 1;
  modconptr->kbd1 = 2;
  modconptr->kbd2 = 2;

  modconptr->yb1 = CONE;
  modconptr->yb2 = CONE;
  modconptr->zb1 = CONE;
  modconptr->zb2 = CONE;

  for (i = 0; i < MAXSIMLOOP; i++)
    {
      structptr->lzcnt[i] = 0;
    }

  structptr->loopzv = 0;

  return (0);
}
