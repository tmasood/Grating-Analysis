/* jairy.f -- translated by f2c (version 19961209).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* DECK JAIRY */
/* Subroutine */ int jairy_(x, rx, c__, ai, dai)
real *x, *rx, *c__, *ai, *dai;
{
    /* Initialized data */

    static integer n1 = 14;
    static real con2 = (float)5.03154716196777;
    static real con3 = (float).380004589867293;
    static real con4 = (float).833333333333333;
    static real con5 = (float).866025403784439;
    static real ak1[14] = { (float).220423090987793,(float)-.1252902427877,(
	    float).0103881163359194,(float)8.22844152006343e-4,(float)
	    -2.34614345891226e-4,(float)1.63824280172116e-5,(float)
	    3.06902589573189e-7,(float)-1.29621999359332e-7,(float)
	    8.22908158823668e-9,(float)1.53963968623298e-11,(float)
	    -3.39165465615682e-11,(float)2.03253257423626e-12,(float)
	    -1.10679546097884e-14,(float)-5.1616949778508e-15 };
    static real ak2[23] = { (float).274366150869598,(float).00539790969736903,
	    (float)-.0015733922062119,(float)4.2742752824875e-4,(float)
	    -1.12124917399925e-4,(float)2.88763171318904e-5,(float)
	    -7.36804225370554e-6,(float)1.87290209741024e-6,(float)
	    -4.75892793962291e-7,(float)1.21130416955909e-7,(float)
	    -3.09245374270614e-8,(float)7.92454705282654e-9,(float)
	    -2.03902447167914e-9,(float)5.26863056595742e-10,(float)
	    -1.36704767639569e-10,(float)3.56141039013708e-11,(float)
	    -9.3138829654843e-12,(float)2.44464450473635e-12,(float)
	    -6.43840261990955e-13,(float)1.70106030559349e-13,(float)
	    -4.50760104503281e-14,(float)1.19774799164811e-14,(float)
	    -3.19077040865066e-15 };
    static real ak3[14] = { (float).280271447340791,(float)
	    -.00178127042844379,(float)4.03422579628999e-5,(float)
	    -1.63249965269003e-6,(float)9.21181482476768e-8,(float)
	    -6.52294330229155e-9,(float)5.47138404576546e-10,(float)
	    -5.2440825180026e-11,(float)5.60477904117209e-12,(float)
	    -6.56375244639313e-13,(float)8.31285761966247e-14,(float)
	    -1.12705134691063e-14,(float)1.62267976598129e-15,(float)
	    -2.46480324312426e-16 };
    static real ajp[19] = { (float).0778952966437581,(float)-.184356363456801,
	    (float).0301412605216174,(float).0305342724277608,(float)
	    -.00495424702513079,(float)-.00172749552563952,(float)
	    2.4313763783919e-4,(float)5.04564777517082e-5,(float)
	    -6.16316582695208e-6,(float)-9.03986745510768e-7,(float)
	    9.70243778355884e-8,(float)1.09639453305205e-8,(float)
	    -1.04716330588766e-9,(float)-9.60359441344646e-11,(float)
	    8.25358789454134e-12,(float)6.36123439018768e-13,(float)
	    -4.96629614116015e-14,(float)-3.29810288929615e-15,(float)
	    2.35798252031104e-16 };
    static real ajn[19] = { (float).0380497887617242,(float)-.245319541845546,
	    (float).165820623702696,(float).0749330045818789,(float)
	    -.0263476288106641,(float)-.00592535597304981,(float)
	    .00144744409589804,(float)2.18311831322215e-4,(float)
	    -4.10662077680304e-5,(float)-4.66874994171766e-6,(float)
	    7.1521880727716e-7,(float)6.52964770854633e-8,(float)
	    -8.44284027565946e-9,(float)-6.44186158976978e-10,(float)
	    7.20802286505285e-11,(float)4.72465431717846e-12,(float)
	    -4.66022632547045e-13,(float)-2.67762710389189e-14,(float)
	    2.36161316570019e-15 };
    static real a[15] = { (float).490275424742791,(float).00157647277946204,(
	    float)-9.66195963140306e-5,(float)1.35916080268815e-7,(float)
	    2.98157342654859e-7,(float)-1.86824767559979e-8,(float)
	    -1.03685737667141e-9,(float)3.28660818434328e-10,(float)
	    -2.5709141063278e-11,(float)-2.32357655300677e-12,(float)
	    9.57523279048255e-13,(float)-1.20340828049719e-13,(float)
	    -2.90907716770715e-15,(float)4.55656454580149e-15,(float)
	    -9.99003874810259e-16 };
    static integer n2 = 23;
    static real b[15] = { (float).278593552803079,(float)-.00352915691882584,(
	    float)-2.31149677384994e-5,(float)4.7131784226356e-6,(float)
	    -1.12415907931333e-7,(float)-2.00100301184339e-8,(float)
	    2.60948075302193e-9,(float)-3.55098136101216e-11,(float)
	    -3.50849978423875e-11,(float)5.83007187954202e-12,(float)
	    -2.04644828753326e-13,(float)-1.10529179476742e-13,(float)
	    2.87724778038775e-14,(float)-2.88205111009939e-15,(float)
	    -3.32656311696166e-16 };
    static integer n1d = 14;
    static integer n2d = 24;
    static integer n3d = 19;
    static integer n4d = 15;
    static integer m1d = 12;
    static integer m2d = 22;
    static integer m3d = 17;
    static integer m4d = 13;
    static real dak1[14] = { (float).204567842307887,(float)
	    -.0661322739905664,(float)-.00849845800989287,(float)
	    .00312183491556289,(float)-2.70016489829432e-4,(float)
	    -6.35636298679387e-6,(float)3.02397712409509e-6,(float)
	    -2.18311195330088e-7,(float)-5.36194289332826e-10,(float)
	    1.1309803562231e-9,(float)-7.43023834629073e-11,(float)
	    4.28804170826891e-13,(float)2.23810925754539e-13,(float)
	    -1.39140135641182e-14 };
    static integer n3 = 19;
    static real dak2[24] = { (float).29333234388323,(float)
	    -.00806196784743112,(float).0024254017233314,(float)
	    -6.82297548850235e-4,(float)1.85786427751181e-4,(float)
	    -4.97457447684059e-5,(float)1.32090681239497e-5,(float)
	    -3.49528240444943e-6,(float)9.24362451078835e-7,(float)
	    -2.44732671521867e-7,(float)6.4930783764891e-8,(float)
	    -1.72717621501538e-8,(float)4.60725763604656e-9,(float)
	    -1.2324905529155e-9,(float)3.30620409488102e-10,(float)
	    -8.89252099772401e-11,(float)2.39773319878298e-11,(float)
	    -6.4801392115345e-12,(float)1.75510132023731e-12,(float)
	    -4.76303829833637e-13,(float)1.2949824110081e-13,(float)
	    -3.5267962221043e-14,(float)9.62005151585923e-15,(float)
	    -2.62786914342292e-15 };
    static real dak3[14] = { (float).284675828811349,(float).0025307307261908,
	    (float)-4.83481130337976e-5,(float)1.84907283946343e-6,(float)
	    -1.01418491178576e-7,(float)7.05925634457153e-9,(float)
	    -5.85325291400382e-10,(float)5.56357688831339e-11,(float)
	    -5.908890947795e-12,(float)6.88574353784436e-13,(float)
	    -8.68588256452194e-14,(float)1.17374762617213e-14,(float)
	    -1.68523146510923e-15,(float)2.55374773097056e-16 };
    static real dajp[19] = { (float).0653219131311457,(float)
	    -.120262933688823,(float).00978010236263823,(float)
	    .0167948429230505,(float)-.00197146140182132,(float)
	    -8.45560295098867e-4,(float)9.42889620701976e-5,(float)
	    2.25827860945475e-5,(float)-2.29067870915987e-6,(float)
	    -3.76343991136919e-7,(float)3.45663933559565e-8,(float)
	    4.29611332003007e-9,(float)-3.58673691214989e-10,(float)
	    -3.57245881361895e-11,(float)2.72696091066336e-12,(float)
	    2.26120653095771e-13,(float)-1.58763205238303e-14,(float)
	    -1.12604374485125e-15,(float)7.31327529515367e-17 };
    static real dajn[19] = { (float).0108594539632967,(float)
	    .0853313194857091,(float)-.315277068113058,(float)
	    -.0878420725294257,(float).0553251906976048,(float)
	    .00941674060503241,(float)-.00332187026018996,(float)
	    -4.11157343156826e-4,(float)1.01297326891346e-4,(float)
	    9.87633682208396e-6,(float)-1.87312969812393e-6,(float)
	    -1.50798500131468e-7,(float)2.32687669525394e-8,(float)
	    1.59599917419225e-9,(float)-2.07665922668385e-10,(float)
	    -1.24103350500302e-11,(float)1.39631765331043e-12,(float)
	    7.3940097115574e-14,(float)-7.328874756275e-15 };
    static real da[15] = { (float).491627321104601,(float).00311164930427489,(
	    float)8.23140762854081e-5,(float)-4.61769776172142e-6,(float)
	    -6.13158880534626e-8,(float)2.8729580465652e-8,(float)
	    -1.81959715372117e-9,(float)-1.44752826642035e-10,(float)
	    4.53724043420422e-11,(float)-3.99655065847223e-12,(float)
	    -3.24089119830323e-13,(float)1.62098952568741e-13,(float)
	    -2.40765247974057e-14,(float)1.69384811284491e-16,(float)
	    8.17900786477396e-16 };
    static real db[15] = { (float)-.277571356944231,(float).0044421283341992,(
	    float)-8.42328522190089e-5,(float)-2.5804031841871e-6,(float)
	    3.42389720217621e-7,(float)-6.24286894709776e-9,(float)
	    -2.36377836844577e-9,(float)3.16991042656673e-10,(float)
	    -4.40995691658191e-12,(float)-5.18674221093575e-12,(float)
	    9.64874015137022e-13,(float)-4.9019057660871e-14,(float)
	    -1.77253430678112e-14,(float)5.55950610442662e-15,(float)
	    -7.1179333757953e-16 };
    static integer n4 = 15;
    static integer m1 = 12;
    static integer m2 = 21;
    static integer m3 = 17;
    static integer m4 = 13;
    static real fpi12 = (float)1.30899693899575;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(), exp(), cos(), sin();

    /* Local variables */
    static real rtrx, temp1, temp2;
    static integer i__, j;
    static real t, e1, e2, f1, f2, ec, cv, tt, ccv, scv;

/* ***BEGIN PROLOGUE  JAIRY */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BESJ and BESY */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (JAIRY-S, DJAIRY-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/*           Daniel, S. L., (SNLA) */
/*           Weston, M. K., (SNLA) */
/* ***DESCRIPTION */

/*                  JAIRY computes the Airy function AI(X) */
/*                   and its derivative DAI(X) for ASYJY */

/*                                   INPUT */

/*         X - Argument, computed by ASYJY, X unrestricted */
/*        RX - RX=SQRT(ABS(X)), computed by ASYJY */
/*         C - C=2.*(ABS(X)**1.5)/3., computed by ASYJY */

/*                                  OUTPUT */

/*        AI - Value of function AI(X) */
/*       DAI - Value of the derivative DAI(X) */

/* ***SEE ALSO  BESJ, BESY */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750101  DATE WRITTEN */
/*   891009  Removed unreferenced variable.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR section.  (WRB) */
/* ***END PROLOGUE  JAIRY */

/* ***FIRST EXECUTABLE STATEMENT  JAIRY */
    if (*x < (float)0.) {
	goto L90;
    }
    if (*c__ > (float)5.) {
	goto L60;
    }
    if (*x > (float)1.2) {
	goto L30;
    }
    t = (*x + *x - (float)1.2) * con4;
    tt = t + t;
    j = n1;
    f1 = ak1[j - 1];
    f2 = (float)0.;
    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + ak1[j - 1];
	f2 = temp1;
/* L10: */
    }
    *ai = t * f1 - f2 + ak1[0];

    j = n1d;
    f1 = dak1[j - 1];
    f2 = (float)0.;
    i__1 = m1d;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + dak1[j - 1];
	f2 = temp1;
/* L20: */
    }
    *dai = -(t * f1 - f2 + dak1[0]);
    return 0;

L30:
    t = (*x + *x - con2) * con3;
    tt = t + t;
    j = n2;
    f1 = ak2[j - 1];
    f2 = (float)0.;
    i__1 = m2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + ak2[j - 1];
	f2 = temp1;
/* L40: */
    }
    rtrx = sqrt(*rx);
    ec = exp(-(*c__));
    *ai = ec * (t * f1 - f2 + ak2[0]) / rtrx;
    j = n2d;
    f1 = dak2[j - 1];
    f2 = (float)0.;
    i__1 = m2d;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + dak2[j - 1];
	f2 = temp1;
/* L50: */
    }
    *dai = -ec * (t * f1 - f2 + dak2[0]) * rtrx;
    return 0;

L60:
    t = (float)10. / *c__ - (float)1.;
    tt = t + t;
    j = n1;
    f1 = ak3[j - 1];
    f2 = (float)0.;
    i__1 = m1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + ak3[j - 1];
	f2 = temp1;
/* L70: */
    }
    rtrx = sqrt(*rx);
    ec = exp(-(*c__));
    *ai = ec * (t * f1 - f2 + ak3[0]) / rtrx;
    j = n1d;
    f1 = dak3[j - 1];
    f2 = (float)0.;
    i__1 = m1d;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	f1 = tt * f1 - f2 + dak3[j - 1];
	f2 = temp1;
/* L80: */
    }
    *dai = -rtrx * ec * (t * f1 - f2 + dak3[0]);
    return 0;

L90:
    if (*c__ > (float)5.) {
	goto L120;
    }
    t = *c__ * (float).4 - (float)1.;
    tt = t + t;
    j = n3;
    f1 = ajp[j - 1];
    e1 = ajn[j - 1];
    f2 = (float)0.;
    e2 = (float)0.;
    i__1 = m3;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	temp2 = e1;
	f1 = tt * f1 - f2 + ajp[j - 1];
	e1 = tt * e1 - e2 + ajn[j - 1];
	f2 = temp1;
	e2 = temp2;
/* L100: */
    }
    *ai = t * e1 - e2 + ajn[0] - *x * (t * f1 - f2 + ajp[0]);
    j = n3d;
    f1 = dajp[j - 1];
    e1 = dajn[j - 1];
    f2 = (float)0.;
    e2 = (float)0.;
    i__1 = m3d;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	temp2 = e1;
	f1 = tt * f1 - f2 + dajp[j - 1];
	e1 = tt * e1 - e2 + dajn[j - 1];
	f2 = temp1;
	e2 = temp2;
/* L110: */
    }
    *dai = *x * *x * (t * f1 - f2 + dajp[0]) + (t * e1 - e2 + dajn[0]);
    return 0;

L120:
    t = (float)10. / *c__ - (float)1.;
    tt = t + t;
    j = n4;
    f1 = a[j - 1];
    e1 = b[j - 1];
    f2 = (float)0.;
    e2 = (float)0.;
    i__1 = m4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	temp2 = e1;
	f1 = tt * f1 - f2 + a[j - 1];
	e1 = tt * e1 - e2 + b[j - 1];
	f2 = temp1;
	e2 = temp2;
/* L130: */
    }
    temp1 = t * f1 - f2 + a[0];
    temp2 = t * e1 - e2 + b[0];
    rtrx = sqrt(*rx);
    cv = *c__ - fpi12;
    ccv = cos(cv);
    scv = sin(cv);
    *ai = (temp1 * ccv - temp2 * scv) / rtrx;
    j = n4d;
    f1 = da[j - 1];
    e1 = db[j - 1];
    f2 = (float)0.;
    e2 = (float)0.;
    i__1 = m4d;
    for (i__ = 1; i__ <= i__1; ++i__) {
	--j;
	temp1 = f1;
	temp2 = e1;
	f1 = tt * f1 - f2 + da[j - 1];
	e1 = tt * e1 - e2 + db[j - 1];
	f2 = temp1;
	e2 = temp2;
/* L140: */
    }
    temp1 = t * f1 - f2 + da[0];
    temp2 = t * e1 - e2 + db[0];
    e1 = ccv * con5 + scv * (float).5;
    e2 = scv * con5 - ccv * (float).5;
    *dai = (temp1 * e1 - temp2 * e2) * rtrx;
    return 0;
} /* jairy_ */

