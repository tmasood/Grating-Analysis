#include <stdio.h>

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
