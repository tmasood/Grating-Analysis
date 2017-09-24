#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double simp_(double (*eval)(double), double a, double b, double eps)
{
  int j, lvl, nrtr[30];
  double da, dx[30], epsp[30], est, est1, est2[30], 
    est3[30], f1, f2[30], f3[30], f4[30], fa, fb, fbp[30], fm, fmp[30], 
    pval[3][30], sum, sx, x2[30], x3[30];
  static int TRUE = 1;

  /****************************************************************
   *
   *     subroutine simp
   *
   *     purpose
   *         definite integration of a real function
   *
   *     usage
   *         f - function to be integrated--user supplied subprogram
   *         a - lower limit of integration
   *         b - upper limit of integration
   *         eps - accuracy desired
   *         value - value of integral
   *
   *     remarks
   *         1.  the values of a,b,eps are unchanged
   *         2.  a valid single precision routine results if the four
   *             real*8 declarations at the beginning of the program
   *             are changed to real
   *         3.  an error stop with an appropriate message occurs if
   *             the algorithm fails to converge
   *         4.  more levels of recursion may be permitted by
   *             changing all occurrences of 30 in the first three
   *             type declarations and in statement number 2 to a
   *             larger integer
   *
   *     subroutines and function subprograms required
   *         none
   *
   *     method
   *         adaptive simpson quadrature using an absolute error
   *         estimate, apportioning the global error equally over all
   *         subintervals, and including a failure to converge stop
   *         reference:  algorithm 182 comm acm 6(1963) 315
   *
   *     october 24, 1968
   ****************************************************************/
  /* the parameter setup for the initial call */
  lvl = 0;
  est = 1;
  da = b - a;
  fa = (*eval)(a);
  fm = 4*(*eval)((a+b)/2);
  fb = (*eval)(b);
  /* recur */
  while( TRUE ){
    lvl += 1;
    dx[lvl-1] = da/3;
    sx = dx[lvl-1]/6;
    f1 = 4*(*eval)(a+dx[lvl-1]/2);
    x2[lvl-1] = a + dx[lvl-1];
    f2[lvl-1] = (*eval)(x2[lvl-1]);
    x3[lvl-1] = x2[lvl-1] + dx[lvl-1];
    f3[lvl-1] = (*eval)(x3[lvl-1]);
    epsp[lvl-1] = eps;
    f4[lvl-1] = 4*(*eval)(x3[lvl-1] + dx[lvl-1]/2);
    fmp[lvl-1] = fm;
    est1 = (fa + f1 + f2[lvl-1])*sx;
    fbp[lvl-1] = fb;
    est2[lvl-1] = (f2[lvl-1] + f3[lvl-1] + fm)*sx;
    est3[lvl-1] = (f3[lvl-1] + f4[lvl-1] + fb)*sx;
    sum = est1 + est2[lvl - 1] + est3[lvl - 1];
    if( fabs( est - sum ) <= epsp[lvl-1]*80 ){
      if( est != 1. ){
	/* done on this level */
	while( TRUE ){
	  lvl -= 1;
	  j = nrtr[lvl-1];
	  pval[j-1][lvl-1] = sum;
	  if( j == 1 )
	    break;
	  if( j == 2 )
	    goto L_100;
	  if( j != 3 )
	    goto L_110;
	  sum = pval[0][lvl-1] + pval[1][lvl-1] + pval[2][lvl-1];
	  if( lvl <= 1 )
	    return(sum);
	}
	nrtr[lvl-1] = 2;
	da = dx[lvl-1];
	fa = f2[lvl-1];
	fm = fmp[lvl-1];
	fb = f3[lvl-1];
	eps = epsp[lvl-1]/3;
	est = est2[lvl-1];
	a = x2[lvl-1];
	continue;
      L_100:
	nrtr[lvl-1] = 3;
	da = dx[lvl-1];
	fa = f3[lvl-1];
	fm = f4[lvl-1];
	fb = fbp[lvl-1];
	eps = epsp[lvl-1]/3;
	est = est3[lvl-1];
	a = x3[lvl-1];
	continue;
      }
    }
    if( lvl >= 30 )
      break;
  L_110:
    nrtr[lvl-1] = 1;
    da = dx[lvl-1];
    fm = f1;
    fb = f2[lvl-1];
    eps = epsp[lvl-1]/3;
    est = est1;
  }
  printf(" error in subroutine simp: failed to converge\n");
  exit( EXIT_FAILURE );
  return(sum);
}

