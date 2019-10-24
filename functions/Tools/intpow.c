#include "typedef.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: intpow.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/*---------------------------------------------------------------------------*\
@
@ int intpow(a,b);
@
@ compute the b-th power of the integer a (int a, int b)
@
\*---------------------------------------------------------------------------*/

int intpow(a,b)
int a,b;
{
int result=1;

  while(b-- > 0) {
    result *= a;
  }
  return(result);
}
