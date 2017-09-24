#include <stdio.h>

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
