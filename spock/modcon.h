#ifndef MODCON_INCLUDE
#define MODCON_INCLUDE

struct MODCON
{
  int kpol;
  float apb1; /* principal branch boundries location  = 0.25 */
  float apb2; /* principal branch boundries location = 0.25 */
  dcomplex vpb1; 
  dcomplex vpb2;
  dcomplex yb1;
  dcomplex yb2;
  dcomplex zb1;
  dcomplex zb2;
  int kbc0;
  int kbc1; /* open boundry eigen conditions = 1 */
  int kbc2; /* open boundry eigen conditions = 1 */
  int kbd1; /* outward only solutions = 2 */
  int kbd2; /* outward only solutions = 2 */
};

#endif
