#ifndef LOCALCOMPLEX_INCLUDE
#define LOCALCOMPLEX_INCLUDE

typedef __complex__ double dcomplex;

extern dcomplex zsqrt(dcomplex);
extern dcomplex zsin( dcomplex);
extern dcomplex zcos( dcomplex);
extern dcomplex zcnjg( dcomplex);
extern double zabsq( dcomplex);
extern double zabs( dcomplex);

extern double log10(double);
extern double pow(double, double);
#endif
