#include <math.h>
#include <complx.h>

/****************************************************
 * The lines with jkb added makes branch cut along
 * the 45 degree line in the third quadrant.
 ****************************************************/

complx zsqrt( complx a )
{
  double dtheta, zm;
  complx v;

  zm = sqrt( sqrt( a.re*a.re+a.im*a.im ) );
  if( a.re == 0.0 ){       /**    The y axis */
    if( a.im > 0.0 ){      /**/
      v.re = zm*M_SQRT1_2; /**    positive y axis */
      v.im = zm*M_SQRT1_2; /**/
    }                      /**/
    else if( a.im < 0.0 ){ /**/
      v.re = zm*M_SQRT1_2; /**    negative y axis */
      v.im = -zm*M_SQRT1_2;/**/
    }                      /**/
    else{                  /**/
      v.re = 0.0;          /**    zero point */
      v.im = 0.0;          /**/
    }
  }
  else if( a.im != 0.0 ){        /* real part not zero im part not zero*/
    dtheta = atan( a.im/a.re );
    if( dtheta <= 0.0 ){        /***/ 
      if( a.re < 0.0 )         /***  second quadrant */
	dtheta += M_PI;       /***/
    }
    else{                        /* dtheta > 0 */
      if( a.re < 0.0 )           /****  third quadrant */
	if( dtheta < 0.25*M_PI )/**jkb**  br cut at 225 degrees */
	  dtheta += M_PI;       /**jkb**/
	else                    /**jkb**/ 
	  dtheta -= M_PI;       /****/
    }
    dtheta *= 0.5;              /* first and forth quadrant */
    v.re = zm*cos( dtheta );  
    v.im = zm*sin( dtheta );
  }
  else if( a.re > 0.0 ){        /* positive x axis */
    v.re = sqrt( a.re );
    v.im = 0.0;
  }
  else{
    v.re = 0.0;               /* negative x axis */
    v.im = sqrt( fabs( a.re ) );
  }
  
  return( v );
}
