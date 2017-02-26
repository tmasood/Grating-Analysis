#include<stdio.h>
#include<math.h>
#include "f2c.h"

extern int cbesj(complex *, real *, integer *,
		 integer *, complex *, integer *, integer *);

extern int cbesjp(complex *, real *, integer *,
		 integer *, complex *, integer *, integer *);

extern int cbesy(complex *, real *, integer *, integer *, complex *,
		 integer *, complex *,  integer *);

extern int cbesyp(complex *, real *, integer *, integer *, complex *,
		 integer *, complex *,  integer *);

extern int cbesi(complex *, real *, integer *, integer *, complex *,
		 integer *, integer *);

extern int cbesip(complex *, real *, integer *, integer *, complex *,
		 integer *, integer *);

extern int cbesk(complex *, real *, integer *, integer *, complex *,
		 integer *, integer *);

extern int cbeskp(complex *, real *, integer *, integer *, complex *,
		 integer *, integer *);

extern int beskp(real*, real*, int*, int*, real*, int*);

int main()
{

  complex z;
  complex cy[10];
  real y[10];
  real z2;
  complex cwrk[10];
  integer kode, n, nz, ierr;
  real fnu;

  z.i = 0.0;
  for (z.r=0.1; z.r<20.1; z.r+=0.1)
    {
      kode = 1;
      fnu = 1.0;
      n = 10;
      cbesj(&z, &fnu, &kode, &n, cy, &nz, &ierr);
      cbesy(&z, &fnu, &kode, &n, cy, &nz, cwrk, &ierr);
      cbesi(&z, &fnu, &kode, &n, cy, &nz, &ierr);
      cbesk(&z, &fnu, &kode, &n, cy, &nz, &ierr);
      cbesjp(&z, &fnu, &kode, &n, cy, &nz, &ierr);
      cbesyp(&z, &fnu, &kode, &n, cy, &nz, cwrk, &ierr);
      cbesip(&z, &fnu, &kode, &n, cy, &nz, &ierr);
      cbeskp(&z, &fnu, &kode, &n, cy, &nz, &ierr);
      printf("complex %f \t %f\n", z.r,cy[0].r );
      z2 = z.r;
      beskp(&z2, &fnu, &kode, &n, y, &nz);
      printf("real    %f \t %f\n", z2,y[0] );
    }
return (0);
}
