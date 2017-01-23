//
//  Copyright 2015, 2016  Taha Masood, Johannes Tausch, Jerome Butler
//  and Gary Evans
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
#include <fstream>
#include <iostream>

#include "ssystem.h"
#include "structure.h"
#include "grating.h"
#include "layer.h"

using namespace std;

int readfile(char* filename, structure *epiptr, grating *gratptr,
	     ssystem *sysptr);
int initialize(structure *epiptr);
int initsheets(int **sheetsleft, int **sheetsright, complex<double> gam,
	       grating *gratptr, structure *epiptr);
complex<double> znewton(structure *epiptr);
int defpnl(structure *epiptr, grating *gratptr);
int defdmn(structure *epiptr, grating *gratptr);
int defpnlst(structure *epiptr, grating *gratptr);
int defdmnst(structure *epiptr, grating *gratptr);
int defpnlslt(structure *epiptr, grating *gratptr);
int defdmnslt(structure *epiptr, grating *gratptr);
int defpnlwg(structure *epiptr, grating *gratptr);
int defdmnwg(structure *epiptr, grating *gratptr);
int indxpnls(grating *gratptr, int *match);
gtpanel *refinepnls(grating *gratptr, int *match);
void gettannrm(grating *gratptr);
int getcollocpts(grating *gratptr);
int calcmoments(structure *epiptr, grating *gratptr, int nrows,
		complex<double> *lambda, complex<double> **at,
		complex<double> **sol0t, complex<double> **sol1t,
		complex<double> **rhst, complex<double> **a1t,
		complex<double> **bt, complex<double> ***mp,
		long **piv);
int subext(complex<double> **mt, complex<double> **extd2nleft,
	    complex<double> **extd2nright, structure *epiptr,
	    grating *gratptr, int *sheetsleft, int *sheetsright);
double infnorm(int n, complex<double> *a);
complex<double> nextgam(int order, complex<double> *z, complex<double> gam,
			int *sheetsleft, int *sheetsright, 
			structure *epiptr, grating *gratptr);
void initcalcp(int qorder);
int dumpnls(grating *gratptr);
complex<double> **getnullspace(int order, int *dimnp, complex<double> *t,
			       double tol);
int printsolution(int order, double x0, double x1, int npts, int npnls, 
		  int nrows, gtpanel *pnls, complex<double> *vleft,
		  complex<double> *vright, complex<double> *extd2nleft,
		  complex<double> *extd2nright, grating *gratptr, 
		  structure *epiptr, int *sheetsleft, int *sheetsright,
		  complex<double> *sol0);

void printsolint(int order, int nptsx, int nptsz, int npnls, int nrows, 
		 gtdomain *dmns, gtpanel *pnls, complex<double> *ns,
		 complex<double> *rhs, complex<double> **sol0t,
		 complex<double> lambda, structure *epiptr,
		 grating *gratptr);

extern "C" {
  // LU factorization
  void  zgetrf_(int *nrow, int *ncol, complex<double> *a, int *lda,
		int *pvt, int *inf);
  // estimate condition number
  void  zgecon_(char *nrm, int *n, complex<double> *a, int *lda,
		double *anorm, double *rcond, complex<double> *wrk,
		double *rwrk, int *inf);
  // LU solve
  void  zgetrs_(char *trans, int *n, int *nrhs, complex<double> *a,
		int *lda, int *pvt, complex<double> *b, int *ldb,
		int *inf);
  // find eigenvalues and left/right eigenvectors of a general
  // complex matrix
  void  zgeev_( char *jobvl, char *jobvr, int *n, complex<double> *a,
		int *lda, complex<double> *w, complex<double> *vl,
		int *ldvl, complex<double> *vr, int *ldvr,
		complex<double> *work, int *lwork, double *rwork,
		int *info);
}

#define ARRAYSIZE 512
#define MAXITER 40

double EPS = 1.0e-10;

int main(int argc, char *argv[])
{
  structure *epiptr;
  grating *gratptr;
  ssystem *sysptr;

  gtpanel *pnls;
  gtdomain *dmns;

  string st;
  string gs;

  char nrm, trans, jobvl, jobvr;

  int ec, loopcnt;
  int status, i;
  int match[ARRAYSIZE];
  int nrows, lda, info;
  int *sheetsleft, *sheetsright;
  int order, ord2, ord4;
  int *pvt;
  int ldvl, ldvr, lwrk;
  int qorder;
  int dimn;
  int j, k;
  int printsolnintvls;
  int npanels;

  complex<double> wgroot; // root of the waveguide without grating
  complex<double> gam;
  complex<double> **mom, *t;
  complex<double> *atmp, *a1tmp, *btmp, *rhstmp, *sol0tmp, *sol1tmp;
  complex<double> *extd2nleft, *extd2nright;
  complex<double> *wrk, *vr, *vl, *eigval;
  complex<double> lambda;
  complex<double> **ns;
  complex<double> gkb;
  complex<double> twoPi;
  complex<double> kapleft, kapright;
  complex<double> I;

  double anorm;
  double rcond;
  double *rwrk;
  double nltol;
  double period;
  double ko;
  double omega, wvl;
  double indx_r, indx_i;
  double loss;
  double astart, aend;

  long *ipv;

  layer *layleftptr, *layrightptr;


  epiptr = new structure;
  gratptr = new grating;
  sysptr = new ssystem;

  atmp = NULL;
  a1tmp = NULL;
  btmp = NULL;
  rhstmp = NULL;
  sol0tmp = NULL;
  sol1tmp = NULL;
  mom = NULL;
  t = NULL;
  extd2nleft = NULL;
  extd2nright = NULL;
  ipv = NULL;
  ns = NULL;
  twoPi = complex<double>(0.0,(2.0 * M_PI));
  I = complex<double>(0.0, 1.0);

  lambda = 1.0;
  qorder = 3;

  astart = -1.0;
  aend = 2.0;
  printsolnintvls = 24000;

  nrm = 'I';
  trans = 'N';
  // Do not compute left and right eigenvectors for mom matrix
  jobvl = 'N';
  jobvr = 'N';

  loopcnt = 0;
  nltol = 1e-10;

  cout << endl;
  cout << endl;

  ec = readfile(argv[argc-1],epiptr,gratptr,sysptr);
  if (ec)
    {
      cout << "Input file not found..." << endl;
      return 0;
    }

  initialize(epiptr);
  cout << endl;
  cout << "     n_eff  " << "               condition "
       << "          info " << endl;
  wgroot = znewton(epiptr);

  cout << " " << endl;
  cout << "Unloaded structure n_eff = " << wgroot << endl;

  st = epiptr->getstructtype();

  if (!st.compare("GRAT"))
    {
      epiptr->gratptr = gratptr;
      gs = epiptr->gratptr->getgratingshape();
      if (!gs.compare("RECT"))
	{
	  initcalcp(qorder);
	  // define panels
	  status = defpnl(epiptr, gratptr);
	  // define domains
	  status = defdmn(epiptr, gratptr);
	}
      else if (!gs.compare("SAWTOOTH"))
	{
	  initcalcp(qorder);
	  // define panels
	  status = defpnlst(epiptr, gratptr);
	  // define domains
	  status = defdmnst(epiptr, gratptr);
	}
      else if (!gs.compare("SLANTED"))
	{
	  initcalcp(qorder);
	  // define panels
	  status = defpnlslt(epiptr, gratptr);
	  // define domains
	  status = defdmnslt(epiptr, gratptr);
	}
      else if (!gs.compare("WG"))
	{
	  initcalcp(qorder);
	  // define panels
	  status = defpnlwg(epiptr, gratptr);
	  // define domains
	  status = defdmnwg(epiptr, gratptr);
	}
      else
	{
	  cout << "grating shape " << gs << " not recognized." << endl;
	  cout << "aborting program ....... " << endl;
	  return 0;
	}

      status = indxpnls(gratptr, match);
      pnls = refinepnls(gratptr, match);
      gettannrm(gratptr);
      nrows = getcollocpts(gratptr);
      // dumpnls(gratptr);

      // get beta guess
      gam = epiptr->getbetaguess();
      epiptr->setbetaguess(gam);
      period = gratptr->getperiod();
      lambda = exp(gam*period);

      order = gratptr->getssspchrmncs();
      ord2 = (2 * order) + 1;
      ord4 = (2 * ord2);
      t = new complex<double>[ord4*ord4];

      eigval = new complex<double>[ord4];
      ldvl = 1;
      ldvr = 1;
      lda = ord4;
      rwrk = new double[2*ord4];
      lwrk = 5*ord4;
      wrk = new complex<double>[lwrk];
      pvt = new int[ord4];
      vr = new complex<double>[ord4];
      vl = new complex<double>[ord4];

      initsheets(&sheetsleft, &sheetsright, gam, gratptr, epiptr);

      while (loopcnt < MAXITER)
	{
	  cout << "gam = " << gam << "\n";
	  calcmoments(epiptr, gratptr, nrows, &lambda, &atmp, &sol0tmp,
		      &sol1tmp, &rhstmp,&a1tmp, &btmp, &mom, &ipv);

	  subext(mom, &extd2nleft, &extd2nright, epiptr, gratptr,
		 sheetsleft, sheetsright);

	  for (i=(ord4*ord4)-1; i>=0; i--)
	    {
	      t[i] = mom[0][i];
	    }

	  // LU factorization
	  anorm = infnorm(ord4, mom[0]);

	  // compute the LU factorization of matrix mom[0]
	  zgetrf_(&lda, &lda, mom[0], &lda, pvt, &info);

	  // estimate the reciprocal of the condition number of matrix mom[0]
	  zgecon_(&nrm, &lda, mom[0], &lda, &anorm, &rcond, wrk, rwrk, &info);

	  cout << "condition = " << rcond << "\n";
	  if (rcond < nltol) break;

	  // solve the system of linear equation with a general matrix mom[0]
	  zgetrs_(&trans, &lda, &lda, mom[0], &lda, pvt, mom[1], &lda, &info);

	  // compute the eigenvalues for matrix mom[1]
	  zgeev_(&jobvl, &jobvr, &lda, mom[1], &lda, eigval, vl, &ldvl, vr,
		 &ldvr, wrk, &lwrk, rwrk, &info);
	  for(i=0; i<ord4; i++)
	    {
	      eigval[i] *= -1;
	    }

	  gam = nextgam(order, eigval, gam, sheetsleft, sheetsright,
			epiptr, gratptr);

	  epiptr->setbetaguess(gam);
	  ko = epiptr->getko();
	  period = gratptr->getperiod();
	  lambda = exp(gam*period);
	  loopcnt++;
	}
      
      if (loopcnt >= MAXITER)
	{
	  cout << endl;
	  cout << " ERROR: program did not converge in " << loopcnt
	       << " iterations" << endl;
	  cout << endl;
	}
      // dump info of calculated eigensolution
      cout << endl;
      cout << "@condition = " << rcond << endl;
      ko = epiptr->getko();
      cout << "Loaded structure n_eff = " << (gam/ko);
      cout << "   loss = " << ((real(gam))*10000.0) << "/cm" << endl;
      cout << endl;
      cout << endl;
    }

  else if (st.compare("WG"))
    {
      cout << endl << endl;
      cout << "Waveguide solution completed ... exiting " << endl;
      cout << endl << endl;
    }

  ns = getnullspace(order, &dimn, t, 100*nltol);

  cout << endl;
  cout << "dimension = " << dimn;
  cout << "  order = " << order << endl;

  layleftptr = epiptr->layerptr;
  layrightptr = epiptr->layerptr;
  while (layrightptr->nextptr != NULL)
    {
      layrightptr = layrightptr->nextptr;
    }

  wvl = epiptr->getwavelength();
  omega = (2 * M_PI)/wvl;
  indx_r = layleftptr->getlayerindex();
  loss = layleftptr->getlayerloss();
  // calculate imaginary part of index from loss
  indx_i = ((loss * wvl)/(4 * M_PI * 10000));
  kapleft = complex<double>(indx_r, indx_i);
  kapleft *= omega;

  indx_r = layrightptr->getlayerindex();
  loss = layrightptr->getlayerloss();
  // calculate imaginary part of index from loss
  indx_i = ((loss * wvl)/(4 * M_PI * 10000));
  kapright = complex<double>(indx_r, indx_i);
  kapright *= omega;

  for (i = 0; i < dimn; i++)
    {
      cout << endl;
      cout << (i+1) << "-th eigensoln, left:" << endl;
      for ( k = -order, j = 0; k <= order; k++, j++)
	{
	  // recalculate transverse propagation factors
	  gkb = sqrt(pow((twoPi*static_cast<double>(k)/period + gam),2.0)
		     + pow((kapleft),2.0));
	  gkb *= static_cast<double>(sheetsleft[order+k]) * I;
	  cout << i <<  "  " << j << "  " <<  k << "  ";
	  cout << abs(ns[i][j]) << "  " <<  ns[i][j] << "  ";
	  cout << gkb << "  " << sheetsleft[order+k] << endl;
	}
          printf("\n %d-th eigensoln, right:\n", i + 1);
	  for (k = -order, j = 2*order+1; k <= order; k++, j++)
	    {
	      //  recalculate transverse propagation factors
	      gkb = sqrt(pow((twoPi*static_cast<double>(k)/period + gam),2.0)
			 + pow((kapright),2.0));
	      gkb *= static_cast<double>(sheetsright[order+k]) * I;
	      cout << i <<  "  " << j << "  " <<  k << "  ";
	      cout << abs(ns[i][j]) << "  " <<  ns[i][j] << "  ";
	      cout << gkb << "  " << sheetsright[order+k] << endl;
	    }
    }

  npanels = gratptr->gettotalpnls();

  printsolution(order, astart, aend, printsolnintvls, npanels, 
		nrows, pnls, &(ns[0][0]), &(ns[0][ord2]), extd2nleft,
		  extd2nright, gratptr, epiptr, sheetsleft,
		  sheetsright, sol0tmp);

  dmns = gratptr->gtdmnptr;

  printsolint(order, 10, 4000, npanels, nrows, dmns, pnls, ns[0],
  	      rhstmp, &sol0tmp, lambda, epiptr, gratptr);

  //  testSol(order, nRows, dmns, nPanels, pnls, NS[0]);


  // delete gratptr;
  // delete epiptr;
  // delete sysptr;

  return 0;
}
