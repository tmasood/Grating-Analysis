#ifndef M_2PI
#define M_2PI 6.28318530717958647692528676656
#endif

#define M_SQRT_2_PI 0.79788456080286535587

/******************************************************************
 *                                                                *
 * Macros : MIN, MAX, ABS, SQR                                    *
 *----------------------------------------------------------------*
 * MIN(A, B) returns the minimum of A and B.                      *
 *                                                                *
 * MAX(A, B) returns the maximum of A and B.                      *
 *                                                                *
 * ABS(A) returns the absolute value of A.                        *
 *                                                                *
 * SQ(A) returns the square of A.                                *
 *                                                                *
 ******************************************************************/
#ifndef MIN
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif

#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#ifndef ABS
#define ABS(A)    ((A < 0) ? -(A) : (A))
#endif

#ifndef SQ
#define SQ(x) ( (x)*(x) )
#endif

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
double adapt(double (*fun)(double , int *, double *), double a, double b, 
             double relerr, double *errest, int *nr, 
	     int *icom, double *rcom);
double adapts(double (*fun)(double , int *, double *), double a, double b, 
             double relerr, double *errest, int *nr, int *ns, 
	     int *icom, double *rcom);
/* Simple version of adapt above */
double Adapt(double (*fun)(double ), double a, double b, 
	     double relerr, double *errest, int *nr );
void itoa(int , char *);
void indexx(int , double *, int *);
double simp_(double (* )(double ), double , double , double );
double seval(int , double , double *, double *, double *, double *, double *);
void spline(int , double *, double *, double *, double *, double *);
double rcseval(double xx, int n, double *x, double *y, double *yp);
void rcspline(int mode, int n, double *x, double *y, double *yp, double *w);
double romb(double (* )(double ), double, double, double );
double powi(double, int );
void floc(double *, double (* )(double), 
	    void (* )(double *, double ,double ));
void reverse(char *);

/* For Airy Function Real */
void airya(double , double *, double *, double *, double *);
void ajyik(double , double *, double *, double *, double *, 
	   double *, double *, double *, double *);

/* Check for allocations */
void IntegerAllocated(int *, char *);
void DoubleAllocated(double *D, char *Dc);

/* Allocate Integers and check if allocated */
int *AllocateInteger(int N, char *Dc);

/* Allocate Doubles and check if allocated */
double *AllocateDouble(int N, char *Dc);

/* Allocate Characters and check if allocated */
char *AllocateCharacter(int N, char *Dc);

void iout_(int , int , int *);
void rout_(int , int , double *);
void charout_(int , int, char *);
