#include <stdio.h>
#include <math.h>
#include  <Jkb.h>

int n;

void zprint(double *array, double z, double e)
{
  double zt;
  zt = 1.0/n;
  ++n;
  printf(" %25.7e %27.5e %15.3e\n",z, zt, e);
  
  return;
}

double fncn(double x)
{
  double y;
  y = sin(M_PI/x);
  return( y );
}

int main()
{
  double Roots, a[5];

  printf(" Root Test\n");
  printf(" y=sin(pi/x) from x=2 to x=1/10.5\n");
  printf("             Computed Zero           True Zero        Predected Error\n");

  a[0] = 2.0;
  a[1] = 1.0/10.5;
  a[2] = 0.1;
  n = 1;
  Roots = Floc(a, fncn, zprint);
  return 0;
}
