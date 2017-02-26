#include<stdio.h>
#include<math.h>
#include "f2c.h"
#include "cmplx.h"

#define SIZE 4
real ccaldet(complex mwave[][SIZE], complex *det)
{
  complex d1, d2, d3, d;
  complex p1, p2, p;

  p1 = c_prod(mwave[2][2], mwave[3][3]);
  p2 = c_prod(mwave[2][3], mwave[3][2]);
  p = c_sub(p1, p2);
  p = c_prod(p, mwave[1][1]);
  d1 = p;

  p1 = c_prod(mwave[2][1], mwave[3][2]);
  p2 = c_prod(mwave[2][2], mwave[3][1]);
  p = c_sub(p1, p2);
  p = c_prod(p, mwave[1][3]);
  d1 = c_add(p, d1);

  p1 = c_prod(mwave[2][2], mwave[3][3]);
  p2 = c_prod(mwave[2][3], mwave[3][2]);
  p = c_sub(p1, p2);
  p = c_prod(p, mwave[1][0]);  
  d2 = p;

  p1 = c_prod(mwave[2][0], mwave[3][2]);
  p2 = c_prod(mwave[2][2], mwave[3][0]);
  p = c_sub(p1, p2);
  p = c_prod(p, mwave[1][3]);
  d2 = c_add(p, d2);
 

  p1 = c_prod(mwave[2][1], mwave[3][3]);
  p2 = c_prod(mwave[2][3], mwave[3][1]);
  p = c_sub(p1, p2);
  p = c_prod(p, mwave[1][0]);  
  d3 = p;


  p1 = c_prod(mwave[2][0], mwave[3][3]);
  p2 = c_prod(mwave[2][3], mwave[3][0]);
  p = c_sub(p1, p2);
  p = c_prod(p, mwave[1][1]);
  d3 = c_sub(d3, p);

  p1 = c_prod(mwave[2][0], mwave[3][1]);
  p2 = c_prod(mwave[2][1], mwave[3][0]);
  p = c_sub(p1, p2);
  p = c_prod(p, mwave[1][3]);
  d3 = c_add(p, d3);

  p1 = c_prod(mwave[0][0], d1);
  p2 = c_prod(mwave[0][1], d2);
  p = c_sub(p1, p2);
  p1 = c_prod(mwave[0][2], d3);
  p = c_add(p, p1);
  d = p;

  *det = d;
  return (0.0);
}
