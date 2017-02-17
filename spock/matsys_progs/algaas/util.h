#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <math.h>

typedef long int integer;
typedef unsigned long int uinteger;
typedef char *address;
typedef short int shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
typedef long int longint;
typedef unsigned long int ulongint;

extern double cabs(complex *);
extern void cdiv(complex *, complex *, complex *);
extern void cexp(complex *, complex *);
extern void csqrt(complex *, complex *);
extern double rimag(complex *);
extern void clog(complex *, complex *);
extern void rcnjg(complex *, complex *);
extern double powdd(doublereal *, doublereal *);
extern double rsign(real*, real*);

#define qbit_clear(a,b) ((a) & ~((ulongint)1 << (b)))
#define qbit_set(a,b)   ((a) |  ((ulongint)1 << (b)))
 
#define TRUE 1
#define FALSE 0

/* max function macro */
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmax(a,b) (doublereal)max(a,b)

/* min function macro */
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)

/* absolute value faction macro */
#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)

#endif
