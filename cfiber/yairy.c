/* yairy.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* DECK YAIRY */
/* Subroutine */ int yairy_(x, rx, c__, bi, dbi)
real *x, *rx, *c__, *bi, *dbi;
{
    /* Initialized data */

    static integer n1 = 20;
    static integer n4d = 14;
    static integer m1d = 19;
    static integer m2d = 18;
    static integer m3d = 17;
    static integer m4d = 12;
    static real fpi12 = (float)1.30899693899575;
    static real spi12 = (float)1.83259571459405;
    static real con1 = (float).666666666666667;
    static real con2 = (float)7.74148278841779;
    static real con3 = (float).364766105490356;
    static integer n2 = 19;
    static real bk1[20] = { (float)2.43202846447449,(float)2.57132009754685,(
	    float)1.02802341258616,(float).341958178205872,(float)
	    .0841978629889284,(float).0193877282587962,(float)
	    .00392687837130335,(float)6.83302689948043e-4,(float)
	    1.14611403991141e-4,(float)1.74195138337086e-5,(float)
	    2.41223620956355e-6,(float)3.24525591983273e-7,(float)
	    4.03509798540183e-8,(float)4.70875059642296e-9,(float)
	    5.35367432585889e-10,(float)5.70606721846334e-11,(float)
	    5.80526363709933e-12,(float)5.76338988616388e-13,(float)
	    5.42103834518071e-14,(float)4.91857330301677e-15 };
    static real bk2[20] = { (float).574830555784088,(float)
	    -.00691648648376891,(float).00197460263052093,(float)
	    -5.24043043868823e-4,(float)1.22965147239661e-4,(float)
	    -2.27059514462173e-5,(float)2.23575555008526e-6,(float)
	    4.15174955023899e-7,(float)-2.84985752198231e-7,(float)
	    8.50187174775435e-8,(float)-1.70400826891326e-8,(float)
	    2.25479746746889e-9,(float)-1.09524166577443e-10,(float)
	    -3.41063845099711e-11,(float)1.11262893886662e-11,(float)
	    -1.75542944241734e-12,(float)1.36298600401767e-13,(float)
	    8.76342105755664e-15,(float)-4.64063099157041e-15,(float)
	    7.7877275873296e-16 };
    static real bk3[20] = { (float).566777053506912,(float).00263672828349579,
	    (float)5.1230335147313e-5,(float)2.10229231564492e-6,(float)
	    1.4221709511389e-7,(float)1.28534295891264e-8,(float)
	    7.28556219407507e-10,(float)-3.45236157301011e-10,(float)
	    -2.11919115912724e-10,(float)-6.56803892922376e-11,(float)
	    -8.14873160315074e-12,(float)3.03177845632183e-12,(float)
	    1.73447220554115e-12,(float)1.67935548701554e-13,(float)
	    -1.49622868806719e-13,(float)-5.15470458953407e-14,(float)
	    8.7574184185783e-15,(float)7.9673555352572e-15,(float)
	    -1.29566137861742e-16,(float)-1.1187879441752e-15 };
    static real bk4[14] = { (float).485444386705114,(float)
	    -.00308525088408463,(float)6.98748404837928e-5,(float)
	    -2.82757234179768e-6,(float)1.59553313064138e-7,(float)
	    -1.12980692144601e-8,(float)9.47671515498754e-10,(float)
	    -9.08301736026423e-11,(float)9.70776206450724e-12,(float)
	    -1.13687527254574e-12,(float)1.43982917533415e-13,(float)
	    -1.95211019558815e-14,(float)2.81056379909357e-15,(float)
	    -4.26916444775176e-16 };
    static real bjp[19] = { (float).134918611457638,(float)-.319314588205813,(
	    float).0522061946276114,(float).0528869112170312,(float)
	    -.0085810075607735,(float)-.00299211002025555,(float)
	    4.21126741969759e-4,(float)8.73931830369273e-5,(float)
	    -1.06749163477533e-5,(float)-1.56575097259349e-6,(float)
	    1.68051151983999e-7,(float)1.89901103638691e-8,(float)
	    -1.81374004961922e-9,(float)-1.66339134593739e-10,(float)
	    1.4295633578081e-11,(float)1.10179811626595e-12,(float)
	    -8.60187724192263e-14,(float)-5.71248177285064e-15,(float)
	    4.08414552853803e-16 };
    static real bjn[19] = { (float).0659041673525697,(float)-.424905910566004,
	    (float).28720974519583,(float).129787771099606,(float)
	    -.0456354317590358,(float)-.010263017598254,(float)
	    .00250704671521101,(float)3.78127183743483e-4,(float)
	    -7.11287583284084e-5,(float)-8.08651210688923e-6,(float)
	    1.23879531273285e-6,(float)1.13096815867279e-7,(float)
	    -1.4623428317631e-8,(float)-1.11576315688077e-9,(float)
	    1.24846618243897e-10,(float)8.18334132555274e-12,(float)
	    -8.07174877048484e-13,(float)-4.63778618766425e-14,(float)
	    4.09043399081631e-15 };
    static real aa[14] = { (float)-.278593552803079,(float).00352915691882584,
	    (float)2.31149677384994e-5,(float)-4.7131784226356e-6,(float)
	    1.12415907931333e-7,(float)2.00100301184339e-8,(float)
	    -2.60948075302193e-9,(float)3.55098136101216e-11,(float)
	    3.50849978423875e-11,(float)-5.83007187954202e-12,(float)
	    2.04644828753326e-13,(float)1.10529179476742e-13,(float)
	    -2.87724778038775e-14,(float)2.88205111009939e-15 };
    static real bb[14] = { (float)-.490275424742791,(float)
	    -.00157647277946204,(float)9.66195963140306e-5,(float)
	    -1.35916080268815e-7,(float)-2.98157342654859e-7,(float)
	    1.86824767559979e-8,(float)1.03685737667141e-9,(float)
	    -3.28660818434328e-10,(float)2.5709141063278e-11,(float)
	    2.32357655300677e-12,(float)-9.57523279048255e-13,(float)
	    1.20340828049719e-13,(float)2.90907716770715e-15,(float)
	    -4.55656454580149e-15 };
    static real dbk1[21] = { (float)2.95926143981893,(float)3.86774568440103,(
	    float)1.80441072356289,(float).578070764125328,(float)
	    .163011468174708,(float).0392044409961855,(float)
	    .00790964210433812,(float).00150640863167338,(float)
	    2.56651976920042e-4,(float)3.93826605867715e-5,(float)
	    5.81097771463818e-6,(float)7.86881233754659e-7,(float)
	    9.93272957325739e-8,(float)1.21424205575107e-8,(float)
	    1.38528332697707e-9,(float)1.50190067586758e-10,(float)
	    1.58271945457594e-11,(float)1.57531847699042e-12,(float)
	    1.50774055398181e-13,(float)1.40594335806564e-14,(float)
	    1.24942698777218e-15 };
    static real dbk2[20] = { (float).549756809432471,(float)
	    .00913556983276901,(float)-.00253635048605507,(float)
	    6.60423795342054e-4,(float)-1.55217243135416e-4,(float)
	    3.00090325448633e-5,(float)-3.76454339467348e-6,(float)
	    -1.33291331611616e-7,(float)2.42587371049013e-7,(float)
	    -8.07861075240228e-8,(float)1.71092818861193e-8,(float)
	    -2.41087357570599e-9,(float)1.53910848162371e-10,(float)
	    2.5646537319063e-11,(float)-9.88581911653212e-12,(float)
	    1.60877986412631e-12,(float)-1.20952524741739e-13,(float)
	    -1.0697827841082e-14,(float)5.02478557067561e-15,(float)
	    -8.68986130935886e-16 };
    static integer n3 = 14;
    static real dbk3[20] = { (float).560598509354302,(float)
	    -.00364870013248135,(float)-5.98147152307417e-5,(float)
	    -2.33611595253625e-6,(float)-1.64571516521436e-7,(float)
	    -2.06333012920569e-8,(float)-4.2774543157311e-9,(float)
	    -1.08494137799276e-9,(float)-2.37207188872763e-10,(float)
	    -2.22132920864966e-11,(float)1.07238008032138e-11,(float)
	    5.71954845245808e-12,(float)7.51102737777835e-13,(float)
	    -3.81912369483793e-13,(float)-1.75870057119257e-13,(float)
	    6.69641694419084e-15,(float)2.26866724792055e-14,(float)
	    2.69898141356743e-15,(float)-2.67133612397359e-15,(float)
	    -6.54121403165269e-16 };
    static real dbk4[14] = { (float).493072999188036,(float)
	    .00438335419803815,(float)-8.37413882246205e-5,(float)
	    3.20268810484632e-6,(float)-1.7566197954827e-7,(float)
	    1.22269906524508e-8,(float)-1.01381314366052e-9,(float)
	    9.63639784237475e-11,(float)-1.02344993379648e-11,(float)
	    1.19264576554355e-12,(float)-1.50443899103287e-13,(float)
	    2.03299052379349e-14,(float)-2.91890652008292e-15,(float)
	    4.42322081975475e-16 };
    static real dbjp[19] = { (float).113140872390745,(float)-.208301511416328,
	    (float).0169396341953138,(float).0290895212478621,(float)
	    -.00341467131311549,(float)-.00146455339197417,(float)
	    1.63313272898517e-4,(float)3.91145328922162e-5,(float)
	    -3.96757190808119e-6,(float)-6.51846913772395e-7,(float)
	    5.9870749526928e-8,(float)7.44108654536549e-9,(float)
	    -6.21241056522632e-10,(float)-6.18768017313526e-11,(float)
	    4.72323484752324e-12,(float)3.91652459802532e-13,(float)
	    -2.74985937845226e-14,(float)-1.9503649776275e-15,(float)
	    1.26669643809444e-16 };
    static real dbjn[19] = { (float)-.018809126006885,(float)-.14779818082614,
	    (float).546075900433171,(float).152146932663116,(float)
	    -.0958260412266886,(float)-.016310273169613,(float)
	    .00575364806680105,(float)7.12145408252655e-4,(float)
	    -1.75452116846724e-4,(float)-1.71063171685128e-5,(float)
	    3.2443558063168e-6,(float)2.61190663932884e-7,(float)
	    -4.03026865912779e-8,(float)-2.76435165853895e-9,(float)
	    3.59687929062312e-10,(float)2.14953308456051e-11,(float)
	    -2.41849311903901e-12,(float)-1.28068004920751e-13,(float)
	    1.26939834401773e-14 };
    static real daa[14] = { (float).277571356944231,(float)-.0044421283341992,
	    (float)8.42328522190089e-5,(float)2.5804031841871e-6,(float)
	    -3.42389720217621e-7,(float)6.24286894709776e-9,(float)
	    2.36377836844577e-9,(float)-3.16991042656673e-10,(float)
	    4.40995691658191e-12,(float)5.18674221093575e-12,(float)
	    -9.64874015137022e-13,(float)4.9019057660871e-14,(float)
	    1.77253430678112e-14,(float)-5.55950610442662e-15 };
    static real dbb[14] = { (float).491627321104601,(float).00311164930427489,
	    (float)8.23140762854081e-5,(float)-4.61769776172142e-6,(float)
	    -6.13158880534626e-8,(float)2.8729580465652e-8,(float)
	    -1.81959715372117e-9,(float)-1.44752826642035e-10,(float)
	    4.53724043420422e-11,(float)-3.99655065847223e-12,(float)
	    -3.24089119830323e-13,(float)1.62098952568741e-13,(float)
	    -2.40765247974057e-14,(float)1.69384811284491e-16 };
    static integer m1 = 18;
    static integer m2 = 17;
    static integer m3 = 12;
    static integer n1d = 21;
    static integer n2d = 20;
    static integer n3d = 19;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(), exp(), cos(), sin();

    /* Local variables */
    static real rtrx, temp1, temp2;
    static integer i__, j;
    static real t, d1, d2, e1, e2, f1, f2, s1, s2, tc, ax, cv, ex, tt;

/* ***BEGIN PROLOGUE  YAIRY */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BESJ and BESY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (YAIRY-S, DYAIRY-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/*           Daniel, S. L., (SNLA) */
/* ***DESCRIPTION */

/*                  YAIRY computes the Airy function BI(X) */
/*                   and its derivative DBI(X) for ASYJY */

/*                                     INPUT */

/*         X  - Argument, computed by ASYJY, X unrestricted */
/*        RX  - RX=SQRT(ABS(X)), computed by ASYJY */
/*         C  - C=2.*(ABS(X)**1.5)/3., computed by ASYJY */

/*                                    OUTPUT */
/*        BI  - Value of function BI(X) */
/*       DBI  - Value of the derivative DBI(X) */

/* ***SEE ALSO  BESJ, BESY */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR section.  (WRB) */
/* ***END PROLOGUE  YAIRY */

/* ***FIRST EXECUTABLE STATEMENT  YAIRY */
    ax = dabs(*x);
    *rx = sqrt(ax);
    *c__ = con1 * ax * *rx;
    if (*x < (float)0.) {
	goto L120;
    }
    if (*c__ > (float)8.) {
	goto L60;
    }
    if (*x > (float)2.5) {
	goto L30;
    }
    t = (*x + *x - (float)2.5) * (float).4;
    tt = t + t;
    j = n1;
    f1 = bk1[j - 1];
    f2 = (float)0.;
    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + bk1[j - 1];
	f2 = temp1;
/* L10: */
    }
    *bi = t * f1 - f2 + bk1[0];
    j = n1d;
    f1 = dbk1[j - 1];
    f2 = (float)0.;
    i__1 = m1d;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + dbk1[j - 1];
	f2 = temp1;
/* L20: */
    }
    *dbi = t * f1 - f2 + dbk1[0];
    return 0;
L30:
    rtrx = sqrt(*rx);
    t = (*x + *x - con2) * con3;
    tt = t + t;
    j = n1;
    f1 = bk2[j - 1];
    f2 = (float)0.;
    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + bk2[j - 1];
	f2 = temp1;
/* L40: */
    }
    *bi = (t * f1 - f2 + bk2[0]) / rtrx;
    ex = exp(*c__);
    *bi *= ex;
    j = n2d;
    f1 = dbk2[j - 1];
    f2 = (float)0.;
    i__1 = m2d;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + dbk2[j - 1];
	f2 = temp1;
/* L50: */
    }
    *dbi = (t * f1 - f2 + dbk2[0]) * rtrx;
    *dbi *= ex;
    return 0;

L60:
    rtrx = sqrt(*rx);
    t = (float)16. / *c__ - (float)1.;
    tt = t + t;
    j = n1;
    f1 = bk3[j - 1];
    f2 = (float)0.;
    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + bk3[j - 1];
	f2 = temp1;
/* L70: */
    }
    s1 = t * f1 - f2 + bk3[0];
    j = n2d;
    f1 = dbk3[j - 1];
    f2 = (float)0.;
    i__1 = m2d;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + dbk3[j - 1];
	f2 = temp1;
/* L80: */
    }
    d1 = t * f1 - f2 + dbk3[0];
    tc = *c__ + *c__;
    ex = exp(*c__);
    if (tc > (float)35.) {
	goto L110;
    }
    t = (float)10. / *c__ - (float)1.;
    tt = t + t;
    j = n3;
    f1 = bk4[j - 1];
    f2 = (float)0.;
    i__1 = m3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + bk4[j - 1];
	f2 = temp1;
/* L90: */
    }
    s2 = t * f1 - f2 + bk4[0];
    *bi = (s1 + exp(-tc) * s2) / rtrx;
    *bi *= ex;
    j = n4d;
    f1 = dbk4[j - 1];
    f2 = (float)0.;
    i__1 = m4d;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + dbk4[j - 1];
	f2 = temp1;
/* L100: */
    }
    d2 = t * f1 - f2 + dbk4[0];
    *dbi = rtrx * (d1 + exp(-tc) * d2);
    *dbi *= ex;
    return 0;
L110:
    *bi = ex * s1 / rtrx;
    *dbi = ex * rtrx * d1;
    return 0;

L120:
    if (*c__ > (float)5.) {
	goto L150;
    }
    t = *c__ * (float).4 - (float)1.;
    tt = t + t;
    j = n2;
    f1 = bjp[j - 1];
    e1 = bjn[j - 1];
    f2 = (float)0.;
    e2 = (float)0.;
    i__1 = m2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	temp2 = e1;
	f1 = tt * f1 - f2 + bjp[j - 1];
	e1 = tt * e1 - e2 + bjn[j - 1];
	f2 = temp1;
	e2 = temp2;
/* L130: */
    }
    *bi = t * e1 - e2 + bjn[0] - ax * (t * f1 - f2 + bjp[0]);
    j = n3d;
    f1 = dbjp[j - 1];
    e1 = dbjn[j - 1];
    f2 = (float)0.;
    e2 = (float)0.;
    i__1 = m3d;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	temp2 = e1;
	f1 = tt * f1 - f2 + dbjp[j - 1];
	e1 = tt * e1 - e2 + dbjn[j - 1];
	f2 = temp1;
	e2 = temp2;
/* L140: */
    }
    *dbi = *x * *x * (t * f1 - f2 + dbjp[0]) + (t * e1 - e2 + dbjn[0]);
    return 0;

L150:
    rtrx = sqrt(*rx);
    t = (float)10. / *c__ - (float)1.;
    tt = t + t;
    j = n3;
    f1 = aa[j - 1];
    e1 = bb[j - 1];
    f2 = (float)0.;
    e2 = (float)0.;
    i__1 = m3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	temp2 = e1;
	f1 = tt * f1 - f2 + aa[j - 1];
	e1 = tt * e1 - e2 + bb[j - 1];
	f2 = temp1;
	e2 = temp2;
/* L160: */
    }
    temp1 = t * f1 - f2 + aa[0];
    temp2 = t * e1 - e2 + bb[0];
    cv = *c__ - fpi12;
    *bi = (temp1 * cos(cv) + temp2 * sin(cv)) / rtrx;
    j = n4d;
    f1 = daa[j - 1];
    e1 = dbb[j - 1];
    f2 = (float)0.;
    e2 = (float)0.;
    i__1 = m4d;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	temp2 = e1;
	f1 = tt * f1 - f2 + daa[j - 1];
	e1 = tt * e1 - e2 + dbb[j - 1];
	f2 = temp1;
	e2 = temp2;
/* L170: */
    }
    temp1 = t * f1 - f2 + daa[0];
    temp2 = t * e1 - e2 + dbb[0];
    cv = *c__ - spi12;
    *dbi = (temp1 * cos(cv) - temp2 * sin(cv)) * rtrx;
    return 0;
} /* yairy_ */

