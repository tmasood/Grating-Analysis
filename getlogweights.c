//
//  Copyright 2015 Quantum Designs LLC, Taha Masood, Johannes Tausch
//  and Jerome Butler
//
//  Permission to use, copy, and distribute this software and its
//  documentation for any purpose with or without fee is hereby granted,
//  provided that the above copyright notice appear in all copies and
//  that both the copyright notice and this permission notice appear
//  in supporting documentation.
//
//  This software is provided "as is" without express or implied warranty
//  to the extent permitted by applicable law.
//

#include <iostream>
#include <complex>
#include <cstdlib>
#include "gtoothpnl.h"

//   getlogweights()
//   quadrature for functions with logarithmic singularity at one endpoint
//
//          int_0^1  f(t) dt  \approx \sum w_i f(t_i)
//
//  integrates functions in span [-ln t, 1, t, .., t^(order-2)] exactly
//  using equispaced nodes. This is not pretty, but it does the job.

void getlogweights( int order, double *t, double *w )
{
  if (order == 2)
    {
      t[0]  = 0.2500000000000000; w[0]  = 0.6483797194839226;
      t[1]  = 0.7500000000000000; w[1]  = 0.3516202805160775;
      return;
    }
  if (order == 3)
    {
      t[0]  = 0.1666666666666667; w[0]  = 0.5220479431787609;
      t[1]  = 0.5000000000000000; w[1]  = -0.0440958863575221;
      t[2]  = 0.8333333333333333; w[2]  = 0.5220479431787612;
      return;
    }
  if (order == 4)
    {
      t[0]  = 0.1250000000000000; w[0]  = 0.4357224019597334;
      t[1]  = 0.3750000000000000; w[1]  = -0.2655005392125327;
      t[2]  = 0.6250000000000000; w[2]  = 0.7238338725458657;
      t[3]  = 0.8750000000000000; w[3]  = 0.1059442647069337;
      return;
    }
  if (order == 5)
    {
      t[0]  = 0.1000000000000000; w[0]  = 0.3862435128159312;
      t[1]  = 0.3000000000000000; w[1]  = -0.5033073845970577;
      t[2]  = 0.5000000000000000; w[2]  = 1.2341277435622533;
      t[3]  = 0.7000000000000000; w[3]  = -0.5033073845970584;
      t[4]  = 0.8999999999999999; w[4]  = 0.3862435128159315;
      return;
    }
  if (order == 6)
    {
      t[0]  = 0.0833333333333333; w[0]  = 0.3459382232179754;
      t[1]  = 0.2500000000000000; w[1]  = -0.6562536160898942;
      t[2]  = 0.4166666666666666; w[2]  = 1.7281322321798149;
      t[3]  = 0.5833333333333333; w[3]  = -1.3312572321798337;
      t[4]  = 0.7499999999999999; w[4]  = 0.8734411160899229;
      t[5]  = 0.9166666666666665; w[5]  = 0.0399992767820147;
      return;
    }
  if (order == 7)
    {
      t[0]  = 0.0714285714285714; w[0]  = 0.3180905851333325;
      t[1]  = 0.2142857142857143; w[1]  = -0.8281615663555619;
      t[2]  = 0.3571428571428571; w[2]  = 2.4914976658889141;
      t[3]  = 0.4999999999999999; w[3]  = -2.9628533693333652;
      t[4]  = 0.6428571428571428; w[4]  = 2.4914976658889061;
      t[5]  = 0.7857142857142856; w[5]  = -0.8281615663555588;
      t[6]  = 0.9285714285714284; w[6]  = 0.3180905851333334;
      return;
    }
  if (order == 8)
    {
      t[0]  = 0.0625000000000000; w[0]  = 0.2939670334329557;
      t[1]  = 0.1875000000000000; w[1]  = -0.9516608097583790;
      t[2]  = 0.3125000000000000; w[2]  = 3.2093134064581408;
      t[3]  = 0.4375000000000000; w[3]  = -4.8759524449573588;
      t[4]  = 0.5625000000000000; w[4]  = 5.0092124069288317;
      t[5]  = 0.6875000000000000; w[5]  = -2.7217855046734711;
      t[6]  = 0.8125000000000000; w[6]  = 1.0253721606187585;
      t[7]  = 0.9375000000000000; w[7]  = 0.0115337519505214;
      return;
    }
  if (order == 9)
    {
      t[0]  = 0.0555555555555556; w[0]  = 0.2756543395565091;
      t[1]  = 0.1666666666666667; w[1]  = -1.0896934218073526;
      t[2]  = 0.2777777777777778; w[2]  = 4.1610112397140959;
      t[3]  = 0.3888888888888889; w[3]  = -7.8722735954933762;
      t[4]  = 0.5000000000000000; w[4]  = 10.0506028760576830;
      t[5]  = 0.6111111111111112; w[5]  = -7.8722735954859102;
      t[6]  = 0.7222222222222223; w[6]  = 4.1610112397061343;
      t[7]  = 0.8333333333333335; w[7]  = -1.0896934218035648;
      t[8]  = 0.9444444444444446; w[8]  = 0.2756543395557818;
      return;
    }
  if (order == 10)
    {
      t[0]  = 0.0500000000000000; w[0]  = 0.2592679930005760;
      t[1]  = 0.1500000000000000; w[1]  = -1.1958972002918902;
      t[2]  = 0.2500000000000000; w[2]  = 5.0607480462821801;
      t[3]  = 0.3500000000000000; w[3]  = -11.2067885369103010;
      t[4]  = 0.4500000000000000; w[4]  = 16.7809532089207170;
      t[5]  = 0.5499999999999999; w[5]  = -16.3328674306877670;
      t[6]  = 0.6499999999999999; w[6]  = 10.8690918895044670;
      t[7]  = 0.7499999999999999; w[7]  = -4.4003435650463585;
      t[8]  = 0.8499999999999999; w[8]  = 1.1693757025432829;
      t[9]  = 0.9499999999999998; w[9]  = -0.0035401073149076;
      return;
    }
  if (order == 11)
    {
      t[0]  = 0.0454545454545455; w[0]  = 0.2461071722313164;
      t[1]  = 0.1363636363636364; w[1]  = -1.3132284678641930;
      t[2]  = 0.2272727272727273; w[2]  = 6.1693656946082180;
      t[3]  = 0.3181818181818182; w[3]  = -15.8516014774454470;
      t[4]  = 0.4090909090909091; w[4]  = 27.7827575048842090;
      t[5]  = 0.5000000000000000; w[5]  = -33.0668008524463770;
      t[6]  = 0.5909090909090909; w[6]  = 27.7827575037278790;
      t[7]  = 0.6818181818181819; w[7]  = -15.8516014761035520;
      t[8]  = 0.7727272727272728; w[8]  = 6.1693656938362773;
      t[9]  = 0.8636363636363638; w[9]  = -1.3132284676298560;
      t[10] = 0.9545454545454547; w[10] = 0.2461071722015227;
      return;
    }
  if (order == 12)
    {
      t[0]  = 0.0416666666666667; w[0]  = 0.2340902332043286;
      t[1]  = 0.1250000000000000; w[1]  = -1.4077529906517861;
      t[2]  = 0.2083333333333333; w[2]  = 7.2255020283155469;
      t[3]  = 0.2916666666666666; w[3]  = -20.8834272533587360;
      t[4]  = 0.3749999999999999; w[4]  = 41.4311499565156040;
      t[5]  = 0.4583333333333333; w[5]  = -57.1682962223100050;
      t[6]  = 0.5416666666666666; w[6]  = 56.7503960212088360;
      t[7]  = 0.6250000000000000; w[7]  = -39.9393445026382990;
      t[8]  = 0.7083333333333334; w[8]  = 19.8018199757045540;
      t[9]  = 0.7916666666666667; w[9]  = -6.3362470477717450;
      t[10] = 0.8750000000000001; w[10] = 1.3045968244900163;
      t[11] = 0.9583333333333335; w[11] = -0.0124870227083048;
      return;
    }
  cout << "order = " << order << "is not supported" << endl;
  exit(1);
}
