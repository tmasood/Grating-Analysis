#include <stdlib.h>
#include <stdio.h>
#include <complx.h>

/*
The matrix f is parititioned with reduced order of 1.
the vector v is extracted from f.

          N   
    |            v     |
    |            v     |
f = |     A      v  B  |
    |            v     | M
    | -----------------|           
    |     C      v  D  |

On input f is MxM. On output f is M-1xM-1. v is Nx1.
*/

int matrixmodify(int M, int N, complx *f, complx *v)
{
  int i, j, c;

  /* Load the vector v */
  if(N >= M){
    printf(" N >= M \n");
    printf("M = %d, N = %d \n",M,N);
    printf("Program Terminated \n");
    exit( EXIT_FAILURE );
  }

  j = 0;
  for(i = M*N ; i < M*N + N; i++){
    *(v+j) = *(f+i);
    j++;
  }
  for(i = N*(M + 1)+1; i < M*(N + 1); i++){
    *(v+j) = *(f + i);
    j++;
  }


  j = 0;

  for(c = 0; c < N ; c++){
    for(i = c*M; i < c*M+N ; i++){
      *(f+j)=*(f+i);
      j++;
    }
    for(i = c*M+N+1; i < (c+1)*M; i++){
      *(f+j)=*(f+i);
      j++;
    }
  } 


  for(c = N+1; c < M ; c++){
    for(i = c*M; i < c*M+N ; i++){
      *(f+j) = *(f+i);
      j++;
    }
    for(i = c*M+N+1; i < (c+1)*M; i++){
      *(f+j) = *(f+i);
      j++;
    }   
  }

  return( j );
}
