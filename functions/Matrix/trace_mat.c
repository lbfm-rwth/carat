#include "typedef.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  trace_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ int trace(mat)
@ matrix_TYP *mat;
@
@ calculates the trace of mat.
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
trace (matrix_TYP *mat)
{
  int i; 
  int erg=0;
  for(i=0; i<mat->cols; i++)
    erg += mat->array.SZ[i][i];
  return(erg);
}



/*{{{}}}*/