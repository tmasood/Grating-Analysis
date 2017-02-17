#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(int argc, char *argv[])
{
  float lambdapl;
  float wvl, ewvl, swvl;
  float effindx, step;
  float h,a1,a2,ep,e1,e2,c;
  float epsilon, ee;
  float ec, term1, term2;

  FILE *fp;

  swvl = atof(argv[1]); 
  ewvl = atof(argv[2]); 
  lambdapl = atof(argv[3]);

  if (swvl <= lambdapl)
    {
      printf("<br><b>error: this material system not valid for \n");
      printf("lambda <= lambdapl\n</b><br>");
      exit(1);
    }

  if (ewvl <= swvl)
    {
      printf("<br><b>Start wavelength should be less than \n");
      printf("End wavelength\n</b><br>");
      exit(1);
    }

  if (ewvl <= (swvl + 0.02))
    {
      printf("<br><b>end wavelength should be 0.03um more than \n");
      printf("start wavelength\n</b><br>");
      exit(1);
    }

  if (ewvl >= 1.8)
    {
      printf("<br><b>end wavelength should be less than \n");
      printf("1.8um \n</b><br>");
      exit(1);
    }

  fp = fopen("/tmp/ingaasp_inp.data","w");
  if(fp == NULL)
    {
      printf("<br><b>Cannot open file </b>\n");
      exit(1);
    }

  step = (ewvl - swvl)/100;
  /* Speed of light 'C', Plank's constant 'H'
     and electron charge 'EC' */
  c = 2.9979e8;
  h = 6.6261e-34;
  ec = 1.6022e-19;
  ep = (h*c)/(lambdapl*1e-6*ec);

  /* All energies are in Electron Volts */
  a1 = 13.3510 - (5.4554 * ep) + (1.2332 * pow(ep,2));
  a2 = 0.7140 - (0.3606 * ep);
          
  e1 = 2.5048;
  e2 = 0.1638;
  wvl = swvl;
  while (wvl < ewvl)
    {
      ee = (h * c)/(wvl * 1e-6 * ec);
      term1 = a1/(1 - (pow((ee/(ep + e1)),2)));
      term2 = a2/(1 - (pow((ee/(ep + e2)),2)));
      epsilon = 1 + term1 + term2;
      effindx = sqrt(epsilon);
      fprintf(fp,"%f \t %f",wvl,effindx);
      fprintf(fp,"\n");
      wvl = wvl + step;
    }
  fclose(fp);  
  return (0);
}
