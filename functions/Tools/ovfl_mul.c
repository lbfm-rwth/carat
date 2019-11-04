#include "typedef.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: ovfl_mul.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/* 
@-------------------------------------------------------------------------
@  int ovfl_mul( a, b);
@ int a,b;
@
@  The routine should catch an integer-overflow.
@  Slow. Unluckily "c" does not provide machine-independent access to the
@  integer hardware overflow-routines.
@  if there is no overflow the function returns the product a*b,
@  otherwise the program exits.
@
@-------------------------------------------------------------------------
 */
int 
ovfl_mul (int a, int b)
{ 
register int result;

  result = a * b;
  if ( result / a != b ) {
    fprintf(stderr,"ovfl_mul: Error: Integer overflow during multiplication\n");
    fprintf(stderr,"left operand: 0x%08x, right operand: 0x%08x, result: 0x%08x\n", a, b, result );
    exit(3);
  }                
  return result;
}

/*{{{}}}*/
