#include <stdio.h>
#include <complx.h>

complx lud_(int nt,int *pivot_index, complx *z)
{
   int i, ij, ik, j, k, kj, kk, m, mj, mk;
   complx t, w;

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
