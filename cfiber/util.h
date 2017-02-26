#ifndef UTIL_H
#define UTIL_H


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

#define qbit_clear(a,b) ((a) & ~((ulongint)1 << (b)))
#define qbit_set(a,b)   ((a) |  ((ulongint)1 << (b)))
 
#define TRUE 1
#define FALSE 0
 
/* table of constant values */
static integer c0 = 0;
static integer c1 = 1;
static integer c2 = 2;
static integer c3 = 3;
static integer c4 = 4;
static integer c5 = 5;
static integer c6 = 6;
static integer c7 = 7;
static integer c8 = 8;
static integer c9 = 9;
static integer c10 = 10;
static integer c11 = 11;
static integer c12 = 12;
static integer c13 = 13;
static integer c14 = 14;
static integer c15 = 15;
static integer c16 = 16;

/* used in cbesj function */
static real hpi = (float)1.57079632679489662;

static complex czero = {(float)0.,(float)0.};

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
