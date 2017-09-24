/* Adapted from "Numerical Recipes" */
void indexx(int n,double *arrin, int *indx)
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
