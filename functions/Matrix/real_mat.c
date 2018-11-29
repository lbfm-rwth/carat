#include "typedef.h"
#include "matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  real_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/*{{{}}}*/
/*{{{  real_mat*/
/*
@-------------------------------------------------------------------------
@ void real_mat( mat, rows, cols);
@ matrix_TYP *mat;
@ int rows, cols;
@
| Erzeugt aus mat eine Matrix mit rows Zeilen und cols Spalten.
| Dabei wird die Matrix passend vergroessert und verkleinert
| und mit Nullen aufgefuellt.
@ Generates form mat a matrix with 'rows' rows and 'cols' columns
@ Therefore the matrix is shrinked or enlarged and filled up with zeros.
@
@-------------------------------------------------------------------------
 */
void real_mat(mat, rows, cols)
matrix_TYP *mat;
int rows, cols;
{
int i,j;
int **Z, **N;

  Z = mat->array.SZ;
  N = mat->array.N;

  if (rows < mat->rows) {
    if( Z ) {
      for (i = rows; i < mat->rows; i++) {
        free( Z[i] );
      }
      Z = (int **)realloc((int *)Z,rows * sizeof(int *));
    }
  
    if( N ) {
      for (i = rows; i < mat->rows; i++) {
        free( N[i] );
      }
      N = (int **)realloc((int *)N,rows * sizeof(int *));
    }
  } else if( rows > mat->rows ) {
    if ( !quick_null_mat( mat ) ) {
      mat->flags.Scalar = FALSE;
    }
    if( Z ) {
      Z = (int **)realloc((int *)Z,rows * sizeof(int *));
      for (i = mat->rows; i < rows; i++) {
        Z[i]= (int  *)calloc(cols , sizeof(int ));
      }
    }
    if( N ) {
      N = (int **)realloc((int *)N,rows * sizeof(int *));
      for (i = mat->rows; i < rows; i++) {
        N[i]= (int  *)calloc(cols , sizeof(int ));
      }
    }
  
  }
  mat->rows = rows;
  mat->array.SZ = Z;
  mat->array.N = N;
  if (mat->cols > cols ) {
    if( Z ) {
      for (i = 0; i < rows;i++) {
        Z[i] = (int *)realloc((int *)Z[i], cols*sizeof(int));
      }
    }
  } else if(mat->cols < cols) {
    if( !quick_null_mat( mat ) ) {
      mat->flags.Scalar = FALSE;
    }
    if( Z ) {
      for (i = 0; i < rows;i++) {
        Z[i] = (int *)realloc((int *)Z[i], cols*sizeof(int));
        for (j = mat->cols; j < cols; j++) {
          Z[i][j] = 0;
        }
      }
    }
    if( N ) {
      for (i = 0; i < mat->rows;i++) {
        N[i] = (int *)realloc((int *)N[i], cols*sizeof(int));
        for (j = mat->cols; j < cols; j++) {
          N[i][j] = 1;
        }
      }
    }
  }
  mat->rows = rows;
  mat->cols = cols;
  mat->flags.Symmetric = mat->flags.Symmetric && (rows == cols);
  mat->flags.Diagonal = mat->flags.Diagonal && (rows == cols);
  mat->flags.Scalar = mat->flags.Scalar && (rows == cols);
  mat->array.SZ = Z;
  mat->array.N = N;
}

/*}}}  */


