#include "f2c.h"
#ifndef CMPLX_INCLUDE
#define CMPLX_INCLUDE

extern complex c_prod(complex, complex);
extern complex cscalar_prod(complex, real);
extern complex c_sub(complex, complex);
extern complex c_add(complex, complex);
extern complex cscalar_div(complex, real);
extern int c_comparege(complex, complex);
extern complex cscalar_div2(real, complex);
extern int c_comparel(complex, complex);
extern int c_compareli(complex, complex);
extern int c_comparel2(complex, complex);
extern int c_compareg(complex, complex);
extern int c_comparegi(complex, complex);
extern void pow_ci(complex *, complex *, integer *);
extern void c_div(complex *, complex *, complex *);
extern void c_sqrt(complex *, complex *);
extern double c_abs(complex *);

#endif
