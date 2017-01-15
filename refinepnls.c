//
//  Copyright 2015  Taha Masood, Johannes Tausch and Jerome Butler
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
#include <cstdlib>
#include "gtoothdmn.h"
#include "grating.h"

using namespace std;

// Refines the orignal panels.
// Returns an array of new panels

gtpanel *refinepnls(grating *gratptr, int *match)
{
  double len, d0, d1;
  double l0, l1, lsq0, lsq1;
  double *x0, *x1;
  double x00, x01, x10, x11;
  double xcoll0, xcoll1;
  double pnlen,dx;

  int i, j, ii, idxa, idxb, k, np, npnls;
  int pnltype;
  int *intbuf2, *intbuf;
  int q;
  int npanels;

  gtpanel *pnlptr, *oldpnl;
  gtpanel *panels, *pnl;
  gtdomain *dmn, *dmnptr;

  // count panels
  pnlptr = gratptr->gtpnlptr;
  dx = gratptr->getdx();
  npnls = 0; // total panels

  for (i = 1, oldpnl = pnlptr; oldpnl != NULL; oldpnl = oldpnl->nextptr, i++)
    {
      x0 = oldpnl->getx0();
      x1 = oldpnl->getx1();
      l0 = (x0[0] - x1[0]);
      l1 = (x0[1] - x1[1]);
      lsq0 = pow(l0,2.0);
      lsq1 = pow(l1,2.0);
      len = (lsq0 + lsq1);
      len = sqrt(len);
      np = static_cast<int>(len/dx) + 1;
      npnls += np;
    }
  // total number of refined panels
  gratptr->settotalpnls(npnls);

  // generate array of refined panels
  npnls = 0;
  pnl = new gtpanel();
  gratptr->gtrefpnlptr = pnl;
  intbuf2 = new int[i];
  intbuf2[0] = 0;

  for (i = 0, oldpnl = pnlptr; oldpnl != NULL; oldpnl = oldpnl->nextptr, i++)
    {
      x0 = oldpnl->getx0();
      x1 = oldpnl->getx1();
      l0 = (x0[0] - x1[0]);
      l1 = (x0[1] - x1[1]);
      lsq0 = pow(l0,2.0);
      lsq1 = pow(l1,2.0);
      len = (lsq0 + lsq1);
      len = sqrt(len);
      // number of refined panels for each panel
      np = static_cast<int>(len/dx) + 1;
      // length of each refined panel
      d0 = (x1[0] - x0[0])/np; // distance in x for refined panel
      d1 = (x1[1] - x0[1])/np; // distance in y for refined panel

      for (ii = 0; ii < np; ii++)
        {
          x00 = x0[0] + ii*d0;
          x01 = x0[1] + ii*d1;
	  pnl->setx0(x00, x01);

          x10 = x0[0] + (ii+1)*d0;
          x11 = x0[1] + (ii+1)*d1;
	  pnl->setx1(x10, x11);
	  
          xcoll0 = 0.5*(x00 + x10);
          xcoll1 = 0.5*(x01 + x11);
	  pnl->setxcoll(xcoll0, xcoll1);

          pnlen = len/np;
	  pnl->setlen(pnlen);

	  pnltype = oldpnl->gettype();
          pnl->settype(pnltype);

	  pnl->nextptr = new gtpanel();
	  pnl = pnl->nextptr;
        }
      npnls += np;
      intbuf2[i+1] = npnls;
    }

  pnl = NULL;

  // index new panels, matching \Gamma_0-\Gamma_2 panels
  // have the same index
  idxa = 0;
  idxb = 0;
  for (i = 0, oldpnl = pnlptr; oldpnl != NULL; oldpnl = oldpnl->nextptr, i++)
    {
      pnltype = oldpnl->gettype();
      if (pnltype != 2)
        { // do nothing for type-2 panels
          if (pnltype > 2)
            {
	      panels = gratptr->gtrefpnlptr;
	      for (q = 0; q < intbuf2[i]; q++)
		{
		  panels = panels->nextptr;
		}
              for (ii = intbuf2[i]; ii < intbuf2[i+1]; ii++)
                {
                  panels->setidxu(idxb);
		  idxb++;
                  panels->setidxv(idxa);
		  idxa++;
		  panels = panels->nextptr;
                }
            }
          else // pnltype < 2
            {
	      panels = gratptr->gtrefpnlptr;
	      for (q = 0; q < intbuf2[i]; q++)
		{
		  panels = panels->nextptr;
		}
              for (ii = intbuf2[i]; ii < intbuf2[i+1]; ii++)
                {
                  panels->setidxu(idxa);
		  idxa++;
                  panels->setidxv(idxa);
		  idxa++;
		  panels = panels->nextptr;
                }

              if (pnltype == 0)
                {
                  k = match[i];
                  idxa -= 2*(intbuf2[i+1] - intbuf2[i]);
		  panels = gratptr->gtrefpnlptr;
		  for (q = 0; q < intbuf2[k]; q++)
		    {
		      panels = panels->nextptr;
		    }
                  for (ii = intbuf2[k]; ii < intbuf2[k+1]; ii++)
                    {
                      panels->setidxu(idxa);
		      idxa++;
                      panels->setidxv(idxa);
		      idxa++;
		      panels = panels->nextptr;
                    }
                }
            }
        }
    }

  // update panel list of domains
  dmnptr = gratptr->gtdmnptr;
  for (dmn = dmnptr; dmn != NULL; dmn = dmn->nextptr)
    {
      npanels = dmn->getnpanels();
      intbuf = new int[npanels];
      for (i = 0; i < npanels; i++)
	{
	  intbuf[i] = dmn->indp[i];
	}

      npnls = 0;
      for (i = 0; i < npanels; i++)
        {
          j = abs(intbuf[i]);
          npnls += intbuf2[j+1] - intbuf2[j];
        }
      dmn->refindp = new int[npnls];
      dmn->refornt = new int[npnls];

      for (k = 0, i = 0; i < npanels; i++)
        {
          j = abs(intbuf[i]);
          np = intbuf2[j+1] - intbuf2[j];
          for (ii = 0 ; ii < np; ii++, k++)
            {
              dmn->refindp[k] = intbuf2[j]+ii;
              dmn->refornt[k] = sgn(intbuf[i]);
            }
        }
      dmn->setnpanels(k);
    } // for dmns

  panels = gratptr->gtrefpnlptr;
  return panels;
}
