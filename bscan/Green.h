void dodiag(int na, int nc, int n, double etatmp, double del, double delc, 
	    double *g, double *h );


double diag(double del, double eta);

void dorev(int na, int nc, int n1, double *g1, double *h1);

void revers(int n, int nc, int nl, int nr, double *a, double *b);

void green(int na, int nc, int n, int *iz, int *ix, int *ir, 
	   double del, double delc, double eta, double *gi, double *hi);

void boupts(int na, int nb, int nc, int n1, int n2, 
	    int *iz1, int *iz2, int *ix1, int *ix2, int *ir1, int *ir2);

complx Gn1(double x_, double z_);
complx Gn2(double x_, double z_);

void GratingFields(int Na, int Nc, int Nb,
		   double Del, double Delc,
		   complx *P1, complx *Q1, complx *P2, complx *Q2,
		   int Nx, int Nz, double *xp, double *zp, 
		   complx g, complx *Fld,
		   double toothEta, double fillEta);

