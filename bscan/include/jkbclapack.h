#include <blaswrap.h>


int zcopy_(int *, complx *, int *, complx *, int *);

int zgecon_(char *, int *, complx *, int *, double *, double *, 
	    complx *, double *, int *);

int zgeru_(int *, int *, complx *, complx *, int *, complx *,
           int *, complx *, int *);

int zgemm_(char *, char *, int *, int *,
           int *, complx *, complx *, int *,
           complx *, int *, complx *, complx *, int *);

int zgemv_(char *, int *, int *, complx *, complx *, int *,
           complx *, int *, complx *, complx *, int *);

int zgeev_(char *, char *, int *,
           complx *, int *, complx *, complx *,
           int *, complx *, int *, complx *,
           int *, double *, int *);

int zgetrf_(int *, int *, complx *,
            int *, int *, int *);
            
int zgetri_(int *, complx *, int *,
            int *, complx *, int *, int *);

int zgesv_(int *, int *, complx *,
           int *, int *, complx *, int *, int *);

int dposv_(char *, int *,int *,double *,int *,double *,int *,int *);
