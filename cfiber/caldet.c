#include<stdio.h>
#include<math.h>
#include "f2c.h"

#define SIZE 4
real caldet(real mwave[][SIZE], real *det)
{
  real d1, d2, d3, d;

  d1 = mwave[1][1]*(mwave[2][2]*mwave[3][3] - mwave[2][3]*mwave[3][2]);
  d1 += mwave[1][3]*(mwave[2][1]*mwave[3][2] - mwave[2][2]*mwave[3][1]);

  d2 = mwave[1][0]*(mwave[2][2]*mwave[3][3] - mwave[2][3]*mwave[3][2]);
  d2 += mwave[1][3]*(mwave[2][0]*mwave[3][2] - mwave[2][2]*mwave[3][0]);

  d3 = mwave[1][0]*(mwave[2][1]*mwave[3][3] - mwave[2][3]*mwave[3][1]);
  d3 -= mwave[1][1]*(mwave[2][0]*mwave[3][3] - mwave[2][3]*mwave[3][0]);
  d3 += mwave[1][3]*(mwave[2][0]*mwave[3][1] - mwave[2][1]*mwave[3][0]);

  d = mwave[0][0]*d1 - mwave[0][1]*d2 + mwave[0][2]*d3;
  *det = d;
  return (0.0);
}
