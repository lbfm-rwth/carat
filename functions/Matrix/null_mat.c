#include "typedef.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: null_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ int null_mat(mat)
@ matrix_TYP *mat;
@ 
@ Checks if the matrix has only 0 as entry.
@ return 1 if yes, 0 otherwise.
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
null_mat (matrix_TYP *mat)
{
   return(save_null_mat(mat));
}
/*{{{}}}*/
/*{{{  save_null_mat*/
/**************************************************************************\
@---------------------------------------------------------------------------
@ int save_null_mat(mat)
@ matrix_TYP *mat;
@ 
@ The same as null_mat
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
save_null_mat (matrix_TYP *mat)
{           
int i,j;

  for ( i=0; i < mat->rows;i++) {
    for ( j=0; j < mat->cols;j++) {
      if ( mat->array.SZ[i][j] != 0 ) return FALSE;
    }
  }
  return TRUE;
}

/*}}}  */
/*{{{  quick_null_mat*/
/**************************************************************************\
@---------------------------------------------------------------------------
@ int quick_null_mat(mat)
@ matrix_TYP *mat;
@ 
@ Checks if the matrix has only 0 as entry using mat->flags.
@ If for example mat->flags.Diagonal = 1 only the diagonal entries
@ are checked to be 0.
@ return 1 if yes, 0 otherwise.
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
quick_null_mat (matrix_TYP *mat)
{   
int i, j;

  if ( mat->flags.Scalar ) {
    return( mat->array.SZ[0][0] == 0 );
  } else if ( mat->flags.Diagonal ) {
    for ( i=0;i < mat->rows;i++) {
      if ( mat->array.SZ[i][i] != 0 ) return FALSE;
    }
    return TRUE;
  } else if ( mat->flags.Symmetric ) {
    for ( i=0;i < mat->rows; i ++ ) {
      for ( j=i; j < mat->cols; j++ ) {
        if ( mat->array.SZ[i][j] != 0 ) {
          return FALSE;
        }
      }
    }
    return TRUE;
  } else {
    return ( save_null_mat(mat) );
  }
}


/*}}}  */
