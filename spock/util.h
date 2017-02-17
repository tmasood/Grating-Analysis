#ifndef UTIL_INCLUDE
#define UTIL_INCLUDE

#include "charmatrix.h"

#define FILESIZE 50
#define MAXCHARS 200
#define MAXKEYWORD 10
#define MAXLAYERS 1000
#define MAXPOINTS 1000
#define MAXLOOP 1000
#define MAXSIMLOOP 20
#define MAXMODES 10

#define ED 2.3025850922940
#define AM 1.0e-14

/* coefficients of power series of sinh function */
#define A0  1.0
#define A2  16.6666666666667e-2
#define A4  83.33333333333e-4
#define A6  198.412698413e-6
#define A8  275.5731922e-8
#define A10 250.52108e-10

#define MAX(A,B)  ( (A) > (B) ? (A) : (B) )
#define MAX3(A,B,C)  ( (C) > MAX(A,B) ? (C) : MAX(A,B) )

extern struct CHARMATRX Smatrix, Tmatrix; /* transfer matrix */

#endif
