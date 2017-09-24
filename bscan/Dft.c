#include <math.h>
#include <complx.h>
#include <jkb.h>

void mfft(int nflt, complx *W )
{
  int i;
  double arg;
  double cfft, sfft;

  W[0] = ftoc(1.0,0.0);
  for( i = 1; i < nflt; i++ ){
    arg = M_2PI*i/nflt;
    sfft = -sin( arg );
    cfft = cos( arg );
    W[i] = ftoc(cfft,sfft);
   }
  return ;
}
/*--------------------------------------------------------*/
void mfftH(int nflt, complx *WH )
{
  int i, N;
  double arg;
  double cfft, sfft;

  N = 2*nflt;
  WH[0] = ftoc(1.0,0.0);
  for( i = 1; i < N; i++ ){
    arg = M_PI*i/nflt;
    sfft = -sin( arg );
    cfft = cos( arg );
    WH[i] = ftoc(cfft,sfft);
   }
  return ;
}
/*--------------------------------------------------------*/
void dft(int nflt, int nfl, complx *W, complx *v, complx *V)
{
  int l, nl, n, n_;
  complx w;

  for( l = 0; l < nflt; l++){
    V[l] = ftoc(0.0,0.0);
    for(n = -nfl; n < nfl; n++){
      n_ = n + nfl;
      nl = abs(n*l) % nflt;
      w = W[nl];
      if(n < 0)
	w = zconj(w);
      V[l] = cadd(V[l],cmul(v[n_],w));
    }
  }
  return;
}
/*--------------------------------------------------------*/
void idft(int nflt, int nfl, complx *W, complx *V, complx *v)
{
  int l, nl, n, n_;
  double N;
  complx w;

  N = (double )nflt;

  for( n = -nfl; n < nfl; n++){
    n_ = n + nfl;
    v[n_] = ftoc(0.0,0.0);
    for( l = 0; l < nflt; l++){
      nl = abs(n*l) % nflt;
      w = W[nl];
      if( n > 0 )
	w = zconj(w);
      v[n_] = cadd(v[n_],cmul(V[l],w));
    }
    v[n_] = cdd(v[n_],N);
  }
  return;
}
/*--------------------------------------------------------*/
void convolve(int nflt, complx *p, complx *h, complx *q)
{
  int l, lp, lmlp;

  for(l = 0; l < nflt; l++){
    p[l].re = 0.0;
    p[l].im = 0.0;
    for(lp = 0; lp < nflt; lp++){
      lmlp = l - lp;
      if(lmlp < 0 )
	lmlp += nflt;
      p[l] = cadd(p[l],cmul(h[lmlp],q[lp]));
    }
  }
  return;
}
/*--------------------------------------------------------*/
void PhaseFactor(int nflt, int nfl, int na, complx *WH, complx *PF)
{
  int n, nn, n_;
  complx w;

  for( n = - nfl; n < nfl; n++) {
    nn = ( abs(n)*(na-1) )%(2*nflt);
    w = WH[nn];
    if( n > 0 )
      w = zconj(w);
    n_ = n + nfl;
    PF[n_] = w;
  }
  return;
}


