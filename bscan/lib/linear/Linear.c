#include <stdio.h>
#include "complx.h"

complex lud_(nt,pivot_index,z)
int nt, pivot_index[];
complex z[];
{
   int i, ij, ik, j, k, kj, kk, m, mj, mk;
   complex t, w;

   w.re = 1.0;
   w.im = 0.0;
   
/* Gaussian elimination with partial pivoting.  */

   pivot_index[nt-1] = 1;

/* Determine the row m containing largest element. Pivot. */

   for (k = 0; k < nt-1; k++){
      m = k;
      for (i = k+1; i < nt; i++){
         ik = i + k*nt;
         mk = m + k*nt;
         if (cabs1(z[ik]) > cabs1(z[mk])) m = i;
         }

/* Check for a nonzero pivot; if all are zero, mat is sing. */

      mk = m + k*nt;
      if (cabs1(z[mk]) == 0.0) return ( ftoc(0.0,0.0) );
      pivot_index[k] = m;
      w = cmul(w,z[mk]);
      
/* Interchange the current row k with the pivot row m. */

      if (m != k){
         pivot_index[nt-1] = -pivot_index[nt-1];
         for (j = k; j < nt; j++){
            mj = m + j*nt;
            kj = k + j*nt;
            t = z[mj];
            z[mj] = z[kj];
            z[kj] = t;
            }
         }

/* Eliminate subdiagonal entries of column k. */

      for (i = k+1; i < nt; i++){
         ik = i + k*nt;
         kk = k*(nt + 1);
         t = cdiv(z[ik],z[kk]);
         z[ik] = cneg(t);
         if (cabs1(t) != 0.0){
            for (j = k+1; j < nt; j++){
               ij = i + j*nt;
               kj = k + j*nt;
               z[ij] = csub(z[ij],cmul(t,z[kj]));
               }
            }
         }
      }
   w = cmul(w,z[nt*nt-1]); 
   
   w.re *= pivot_index[nt-1];
   w.im *= pivot_index[nt-1];

   return( w );
}   
/*********************************************************/
void iout_(int nrow, int ncol, int *mat)
   {
   int i, j, ij;

   for( i = 0; i < nrow; i++){
      for(j = 0; j < ncol; j++){
         ij = i + j*nrow;
         printf(" %3d ",*(mat+ij));
         }
      printf("\n");
      }
   return;   
   }         
/*********************************************************/
void rout_(int m, int n, double *a)
   {
   int i,j,k,ik,js,je;

   for( j = 0; j< n ; j += 10){
      printf(" %5d\n",j);
      js = j;
      je = j + 10;
      if(je >= n)
         je = n;
      for( i = 0; i < m; i++){
         for(k = js ; k < je ; k++){
            ik = i + k*m;
            printf(" %9.5f",*(a+ik));
            }
         printf("\n");   
         }
      }
   return;
   }
/*********************************************************/
void cout_(m, n, a)
int m, n;
complex a[];
{
   int i,j,k,ik,js,je;

   for( j = 0; j< n ; j += 5){
      printf(" %5d\n",j);
      js = j;
      je = j + 5;
      if(je >= n)
         je = n;
      for( i = 0; i < m; i++){
         for(k = js ; k < je ; k++){
            ik = i + k*m;
            printf(" %10.3e %10.3e ",a[ik].re,a[ik].im);
            }
         printf("\n");   
         }
      }
   return;
}
/*********************************************************/
/* Adapted from "Numerical Recipes" */
void indexx(n,arrin,indx)
int n,indx[];
double arrin[];
{
   int l,j,ir,indxt,i;
   double q;

   for ( j = 0 ; j < n ; j++ ) 
      indx[j]=j;
   l=n >> 1;
   ir=n-1;
   for (;;) {
      if (l > 0){
         --l;
         indxt=indx[l];
         q=arrin[indxt];
         }
      else {
         indxt=indx[ir];
         q=arrin[indxt];
         indx[ir]=indx[0];
         --ir;
         if (ir == 0) {
            indx[0]=indxt;
            return;
            }
         }
      i=l;
      j=(l << 1) + 1;
      while (j <= ir) {
         if (j < ir )
            if( arrin[indx[j]] < arrin[indx[j+1]]) 
               j++;
         if (q < arrin[indx[j]]){
            indx[i]=indx[j];
            i = j;
            j = j << 1;
            }
         else 
            j=ir+1;
         }
      indx[i]=indxt;
      }
}
