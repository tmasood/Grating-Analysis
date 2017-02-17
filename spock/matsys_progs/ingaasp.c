#include<stdio.h>
#include<math.h>

int main()
{
      long double xlambda;
      long double effindx;
      long double h,a1,a2,ep,e1,e2,c;
      long double epsilon,e;
      long double ec, term1, term2;

      xlambda = 0.93L;
      while (xlambda < 1.8L)
	{
	  /* All energies are in Electron Volts */
	  ep = 1.32L;
	  a1 = 13.3510L - (5.4554L * ep) + (1.2332L * powl(ep,2));
	  a2 = 0.7140L - (0.3606L * ep);
	  
	  e1 = 2.5048L;
	  e2 = 0.1638L;

	  /* Speed of light 'C', Plank's constant 'H'
	     and electron charge 'EC' */
	  c = 2.9979e8;
	  h = 6.6261e-34;
	  ec = 1.6022e-19;
	  e = (h*c)/(xlambda*1e-6*ec);

	  term1 = a1/(1 - (powl((e/(ep + e1)),2)));
	  term2 = a2/(1 - (powl((e/(ep + e2)),2)));
	  epsilon = 1 + term1 + term2;
	  effindx = sqrtl(epsilon);
	  printf("%Lg \t %Lg \n",xlambda, effindx);
	  xlambda = xlambda + 0.001L;
	}
}

