#ifndef M_2PI
#define M_2PI 6.28318530717958647692528676656
#endif

#define SQ(x) ( (x)*(x) )

/* Matrix and vector devitions */
typedef struct{
  double A;  double B;
  double C;  double D;
} Matrix;

typedef struct{
  double a;
  double b;
} Vector;
/********************************************************************/
void itoa(int , char *);
void indexx(int , double *, int *);
double simp_(double (* )(double ), double , double , double );
double seval(int , double , double *, double *, double *, double *, double *);
void spline(int , double *, double *, double *, double *, double *);
double romb(double (* )(double ), double, double, double );
double powi(double, int );
double floc(double *, double (* )(double), void (* )(double *, double ,double ));
void reverse(char *);

/* For Airy Function Real */
void airya(double , double *, double *, double *, double *);
void ajyik(double , double *, double *, double *, double *, 
	   double *, double *, double *, double *);
