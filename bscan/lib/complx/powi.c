#define IS_ODD(j)	((j) & 1 )

double powi( double x, int n )
{
  double p;

  if( x == 0 )
    return( 0. );

  if( n < 0 ){
    n = -n;
    x = 1/x;
  }

  p = IS_ODD(n) ? x : 1;

  while( n >>= 1 ){
    x *= x;
    if( IS_ODD(n) )
      p *= x;
  }
  return(p);
}
