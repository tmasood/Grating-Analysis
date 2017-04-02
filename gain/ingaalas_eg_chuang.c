#include<stdio.h>
#include <stdlib.h>
#include<math.h>

int main(int argc, char *argv[])
{
  double x, y, Eg, wavelength;
  double eps, epsxx, epsyy, epszz, a0, a;
  double aInAs, aGaAs, aAlAs;
  double C11, C12;
  double C11InAs, C11GaAs, C11AlAs;
  double C12InAs, C12GaAs, C12AlAs;
  double delEc, Peps, Qeps;
  double ac, acInAs, acGaAs, acAlAs;
  double av, b, avInAs, avGaAs, avAlAs;
  double bInAs, bGaAs, bAlAs;
  double EvInAs, EvGaAs, EvAlAs, Evh, Ev;
  double EcInAs, EcGaAs, EcAlAs, Ech, Ec;
  double delEhh, delElh, Echh, Eclh;
  double Egstrain, Ecb, Ecw, Evb, Evw;
  double condbandoffset;

  /*  x = atof(argv[1]);
      y = atof(argv[2]); */
  /* x = 0.1247; Ga QW */
  /* y = 0.1653; Al QW */ 
  a0 = 5.8688; /* lattice constant of InP */
  aInAs = 6.0584; /* lattice constant of InAs */
  aGaAs = 5.6533; /* lattice constant of GaAs */
  aAlAs = 5.66; /* lattice constant of AlAs */

  C11InAs = 8.329;
  C11GaAs = 11.879;
  C11AlAs = 12.5;

  C12InAs = 4.526;
  C12GaAs = 5.376;
  C12AlAs = 5.34;  

  acInAs = -5.08;
  acGaAs = -7.17;
  acAlAs = -5.64;  

  avInAs = 1.00;
  avGaAs = 1.16;
  avAlAs = 2.47;  

  bInAs = -1.8;
  bGaAs = -1.7;
  bAlAs = -1.5;

  EvInAs = 0.441;
  EvGaAs = 0.111;
  EvAlAs = -0.4245;

  EcInAs = 0.801;
  EcGaAs = 1.53;
  EcAlAs = 2.5255;

  for (x=0; x <=1; x=x+0.1)
    {
      for (y=0; y<=1; y=y+0.1)
        {
	  a = (aInAs*(1-x-y)) + (aGaAs*x) + (aAlAs*y);

	  /* C11, C12 - elastic stiffness constants */
	  C11 = (C11InAs*(1-x-y)) + (C11GaAs*x) + (C11AlAs*y);
	  C12 = (C12InAs*(1-x-y)) + (C12GaAs*x) + (C12AlAs*y);

	  ac = (acInAs*(1-x-y)) + (acGaAs*x) + (acAlAs*y);
	  av = (avInAs*(1-x-y)) + (avGaAs*x) + (avAlAs*y);

	  /* valence-band  shear deformation potential */
	  b = (bInAs*(1-x-y)) + (bGaAs*x) + (bAlAs*y);

	  /* valence-band position */
	  Evh = (EvInAs*(1-x-y)) + (EvGaAs*x) + (EvAlAs*y);

	  /* conduction-band position */
	  Ech = (EcInAs*(1-x-y)) + (EcGaAs*x) + (EcAlAs*y);

	  Eg = 0.36 + 2.093*y + 0.629*x + 0.577*y*y + 0.436*x*x
	    + 1.013*x*y - 2.0*x*y*(1 - x - y);
	  wavelength = 1.24/Eg;
	  /* printf ("Eg =  %f - %f\n", Eg, wavelength); */

	  /* The strain in the plain of the epitaxial growth */
	  /* a - lattice constant of quaternary epitaxial layer */
	  /* a0 - lattice constant of the substrate InP */
	  eps = (a0 - a)/a;
	  epsxx = eps;
	  epsyy = eps;

	  /* The strain in the perpendicular direction */
	  /* C12 and C11 are the elastic stiffness constants */
	  epszz = -2*(C12/C11)*eps;

	  /* The conduction band is shifted by the energy */
	  /* ac - conduction band hydrostatic deformation potential */
	  /* av - valence band hydrostatic deformation potential */
	  /* b - valence band shear deformation potential */
	  delEc = ac*(epsxx + epsyy + epszz);

	  Peps = -av*(epsxx + epsyy + epszz);
	  Qeps = -(b/2)*(epsxx + epsyy - 2*epszz);

	  /* valence bands are shifted by */
	  delEhh = -Peps - Qeps;
	  delElh = -Peps + Qeps;

	  /* strained bandgaps */
	  Echh = Eg + delEc - delEhh; 
	  Eclh = Eg + delEc - delElh;

	  if (epsxx < 0.0)
	    {
	      Egstrain = Echh;
	    }
	  else
	    {
	      Egstrain = Eclh;
	    }
	  wavelength = 1.24/Egstrain;
	  if ((x+y)>=0.9)
	    {
	      break;
	    }

	  if (epsxx < 0.0)
	    {
	      Ev = Evh + delEhh;
	    }
	  else
	    {
	      Ev = Evh + delElh;
	    }

	  Ec = Ech + delEc;
	  Ecb = 1.35; /* InP conduction band position */
	  Evb = 0.00; /* InP valence band position */
	  Ecw = Ec;
	  Evw = Ev;

	  condbandoffset = (Ecb - Ecw) / ((Evw - Evb) + (Ecb - Ecw));
	  printf ("%f \t %f \t %f \t %f \t %f \n",x,y,epsxx,Egstrain,condbandoffset);
	}
    }

  return 0;
}
