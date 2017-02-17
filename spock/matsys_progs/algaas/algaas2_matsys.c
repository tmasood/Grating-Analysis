#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

int main(int argc, char *argv[])
{
  real p, wvl;
  real a, b1, b2, b11, b22;
  real c, d, e1, ec, eic, einf, e1c;
  real t, g, f, fi;
  real e0, ei, d0, d1;
  real ek, ee, eei, ep2e0, ep1e0;
  real ep1e01, ep1e02;
  complex arg1, arg2;
  complex q1, q2, q3;
  complex ce00, cec, ce, ce0, cexp0, cexp00;
  complex ctmp, ctmp1, ctmp2;
  
  real ep1e11, ep1e12, ep1e1, ep2e11, ep2e12;
  real ep2e1, ep2ec, ep1ec, eq1, ep2ei, epsilon1;
  real epsilon2, neff;

  wvl = atof(argv[1]); 
  p = atof(argv[2]);

  /*parameters calculated using equation  */
  /*Qi(x)=Q0+(Q1*(P-.5))+(Q2*(P-.5)^2) */
  a=4.264374+(4.402701*(p-0.5))+(9.160232*((p-0.5)*(p-0.5)));
  b1=8.256268+(-0.585350*(p-0.5))+(-1.901850*((p-0.5)*(p-0.5)));
  b2=0.462866+(0.012206*(p-0.5))+(1.047697*((p-0.5)*(p-0.5)));
  b11=19.788257+(-3.706511*(p-0.5))+(2.990894*((p-0.5)*(p-0.5)));
  b22=0.539077+(-0.172307*(p-0.5))+(1.031416*((p-0.5)*(p-0.5)));
  c=3.078636+(-1.598544*(p-0.5))+(-0.742071*((p-0.5)*(p-0.5)));
  d=29.803617+(-22.124036*(p-0.5))+(-57.863105*((p-0.5)*(p-0.5)));
  e1=3.212414+(0.804397*(p-0.5))+(0.228317*((p-0.5)*(p-0.5)));
  ec=4.724383+(0.024499*(p-0.5))+(0.030653*((p-0.5)*(p-0.5)));
  eic=3.160138+(0.138046*(p-0.5))+(1.066214*((p-0.5)*(p-0.5)));
  einf=-0.494941+(-0.030802*(p-0.5))+(0.486201*((p-0.5)*(p-0.5)));
  e1c=6.413109+(0.571119*(p-0.5))+(-0.735610*((p-0.5)*(p-0.5)));
  t=0.263476+(0.090532*(p-0.5))+(-0.254099*((p-0.5)*(p-0.5)));
  g=0.147517+(-0.068764*(p-0.5))+(0.047345*((p-0.5)*(p-0.5)));
  f=1.628226+(-1.682422*(p-0.5))+(-2.081273*((p-0.5)*(p-0.5)));
  fi=0.507707+(-0.070165*(p-0.5))+(-0.122169*((p-0.5)*(p-0.5)));
  
  /* {parameters calculated using equation */
  /* Qi(x)=Q0+(Q1*P)+(Q2*P^2) */
  e0=1.425000+(1.155000*p)+(0.370000*(p*p));
  ei=1.734000+(0.574000*p)+(0.055000*(p*p));
  d0=0.340000+(0.0*p)+(0.0*(p*p));
  d1=0.230000+(-0.030000*p)+(0.0*(p*p));
  ek=(4.135701327e-15 * 2.99792458e14)/wvl;
  ee=ek * ek;
  eei=1.0/ee;
  /* real part of dielectric function: direct edge (e0) */
  arg1.r = 1.0 + (ek/e0);
  arg1.i = 0.0;
  arg2.r = 1.0 - (ek/e0);
  arg2.i = 0.0;
  csqrt(&q1, &arg2);
  csqrt(&q2, &arg1);

  ctmp1.r = (2.0 - q2.r - q1.r);
  ctmp1.i = (0.0 - q2.i - q1.i);
  ctmp2.r = sqrt(e0);
  ctmp2.i = 0.0;
  ctmp.r = ((ctmp2.r * ctmp1.r) - (ctmp2.i * ctmp1.i));
  ctmp.i = ((ctmp2.r * ctmp1.i) + (ctmp2.i * ctmp1.r));
  ep1e01 = a * eei * ctmp.r;

  arg1.r = 1.0 + (ek/(e0 + d0));
  arg1.i = 0.0;
  arg2.r = 1.0 - (ek/(e0 + d0));
  arg2.i = 0.0;
  /* ep1e02=a*eei*(0.5*sqrt(e0+d0)*(2.0-sqrt(arg1)-sqrt(arg2))); */
  csqrt(&q1, &arg2);
  csqrt(&q2, &arg1);
  ep1e02 = a * eei * 0.5*(sqrt(e0+d0))*(2.0 - q2.r - q1.r);
  ep1e0 = ep1e01 + ep1e02;

  /* imaginary part of dielectric function: direct edge (e0) */
  /* taking the real part is like a heavyside function */
  /* ep2e0 = ((2*sqrt(abs(arg1))+sqrt(abs(arg2)))*eei)*a/3.0; */
  arg1.r = ek - e0;
  arg1.i = 0.0;
  arg2.r = ek - e0 - d0;
  arg2.i = 0.0;
  csqrt(&q1, &arg1);
  q3.r = (2 * q1.r);
  csqrt(&q1, &arg2);
  ep2e0 = (q3.r + q1.r) * eei * a / 3.0;

  /* real part of the dielectric function: l [k=pi/a] (e1) */
  /* add a linewidth gam */
  /* cec=1.0/(1.0-((ce00/e1c)*(ce00/e1c))); */
  ce00.r = ek;
  ce00.i = t;
  q1.r = ce00.r/e1c;
  q1.i = ce00.i/e1c;
  q2.r = (q1.r * q1.r) - (q1.i * q1.i);
  q2.i = (q1.r * q1.i) + (q1.i * q1.r);
  q3.r = 1.0;
  q3.i = 0.0;
  q1.r = q3.r - q2.r;
  q1.i = q3.i - q2.i;
  cec.r = 1.0/q1.r;
  cec.i = 1.0/q1.i;

  ce.r = (ce00.r * ce00.r) - (ce00.i * ce00.i);
  ce.i = (ce00.r * ce00.i) + (ce00.i * ce00.r);

  q1.r = ((e1 + d1)*(e1 + d1));
  q1.i = 0.0;
  cdiv(&ce0, &ce, &q1);

  q1.r = (e1 * e1);
  q1.i = 0.0;
  cdiv(&q2, &ce, &q1);
  ce = q2;

  q1.r = min(exp(-f*(ek-e1)),1.0);
  q1.i = 0.0;
  cdiv(&cexp0, &q1, &ce);

  q1.r = min(exp(-f*(ek-e1-d1)),1.0);
  q1.i = 0.0;
  cdiv(&cexp00, &q1, &ce0);

  q1.r = 1.0 - ce.r;
  q1.i = 0.0 - ce.i;
  q2.r = (q1.r * cec.r) - (q1.i * cec.i);
  q2.i = (q1.r * cec.i) + (q1.i * cec.r);
  clog(&q1, &q2);
  q2.r = (q1.r * cexp0.r) - (q1.i * cexp0.i);
  q2.i = (q1.r * cexp0.i) + (q1.i * cexp0.r);
  ep1e11 = -b1 * q2.r;

  q1.r = 1.0 - ce0.r;
  q1.i = 0.0 - ce0.i;
  q2.r = (q1.r * cec.r) - (q1.i * cec.i);
  q2.i = (q1.r * cec.i) + (q1.i * cec.r);
  clog(&q1, &q2);
  q2.r = (q1.r * cexp00.r) - (q1.i * cexp00.i);
  q2.i = (q1.r * cexp00.i) + (q1.i * cexp00.r);
  ep1e12 = -b2 * q2.r;

  ep1e1 = ep1e11 + ep1e12;

  /* imaginary part of dielectric function: l [k=pi/a] (e1) */
  /*  don't allow a negative imaginary part */
  arg2.r = e1 + d1 - ce00.r;
  arg2.i = 0.0;
  csqrt(&q1, &arg1);
  q2.r = cexp0.r * (b1 - (q1.r * b11));
  q2.i = 0.0;
  q3.r = max(q2.r, 0.0);
  q3.i = 0.0;
  ep2e11 = q3.r * M_PI;

  csqrt(&q1, &arg2);
  q2.r = cexp00.r * (b2 - (q1.r * b22));
  q2.i = 0.0;
  q3.r = max(q2.r, 0.0);
  q3.i = 0.0;
  ep2e12 = q3.r * M_PI;
  ep2e1 = max((ep2e11 + ep2e12),0.0);
        
  /* real part of dielectric function: (ec) */
  ce.r = 1.0 - (ee/(ec*ec));
  ce.i = 0.0;
  ep2ec = c/((ce.r * ce.r) + (ee * ((g/ec) * (g/ec))));
  ep1ec = (ce.r * ep2ec);

  /* imaginary part of dielectric function: (ec) */
  ep2ec = (ek * g/ec) * ep2ec;
        
  /* imaginary part of dielectric function: */
  /* indirect edge k=[2*pi/a] (x->gamma) */
  eq1 = ek - ei;
  ep2ei = (d * eei *(eq1 * eq1)) * min(exp(-fi*(ek-eic)),1.0);

  /* do all the heavyside functions */
  if (ek < e0)
    {
      ep2e0 = 0.0;
      ep2e1 = 0.0;
      ep2ec = 0.0;
      ep2ei = 0.0;
    }
  if (ek < e1)
    {
      ep2e1=0.0;
    }
        
  /* multiply the imaginary part by heavyside function. */
  /* Index is the real part of the square root */
  /* of the complex dielectric function */
  epsilon1 = einf + ep1e0 + ep1e1 + ep1ec;
  epsilon2= ep2e0 + ep2e1 + ep2ec + ep2ei;
  arg1.r = epsilon1;
  arg1.i = epsilon2;
  csqrt(&q1, &arg1);
  neff = q1.r;
  printf ("<th align=right> %f </th> \n", neff);
  return (0);
}
