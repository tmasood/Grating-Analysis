#include <complx.h>
#include <stdio.h>

void cout_(int m, int n, complx *a)
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
