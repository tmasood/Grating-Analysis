#include <complx.h>


void vectormodify(int M, int N, complx *v)
{
/*****************************************
 * This algorithm shifs and changes the sign
 * of the solution vector v.  
 *
 * M is the vector length.
 * N elements are shifted down.
 * The Nth elements is unity.
 *
 * M-N elements are shifted down with a
 * sign change.
 * The Nth elements is (1.0,0.0).
 * The top N elements reversed in sign.
 *
 ***************************************/ 
  int i;

  for(i = M-1; i > N ; i--){
    (*(v+i)).re = -(*(v+i-1)).re;
    (*(v+i)).im = -(*(v+i-1)).im;
  }

  (*(v+N)).re = 1.0;
  (*(v+N)).im = 0.0;

  for(i = N-1; i >= 0; i--){
    (*(v+i)).re = -(*(v+i)).re;
    (*(v+i)).im = -(*(v+i)).im;
  }   

}
