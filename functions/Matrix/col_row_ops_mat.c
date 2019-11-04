#include "typedef.h"
#include "matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: col_row_ops_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/*{{{}}}*/
/*{{{  row_per*/
/*
@---------------------------------------------------------------------
@   void row_per( M, i, j );
@   matrix_TYP *M;
@   int i , j;
@
@ swaps the i. and j. row by swapping the pointers of the 2-dim. array
@---------------------------------------------------------------------
 */
void 
row_per (matrix_TYP *M, int i, int j)
{
int *merk;

  merk = M->array.SZ[i];
  M->array.SZ[i] = M->array.SZ[j];
  M->array.SZ[j] = merk;
  if ( M->array.N != NULL )
  {
    merk = M->array.N[i];
    M->array.N[i] = M->array.N[j];
    M->array.N[j] = merk;
  }
  M->flags.Symmetric = FALSE;
}

/*}}}  */
/*{{{  col_per*/
/*
@------------------------------------------------------------------------
@ void col_per(M, i, j);
@ 
@ same as row_per() for columns.
@------------------------------------------------------------------------
 */
void 
col_per (matrix_TYP *M, int i, int j)
{
int k;
int merk;

  for(k=0; k < M->rows; k++)
  {
     merk = M->array.SZ[k][i];
     M->array.SZ[k][i] = M->array.SZ[k][j];
     M->array.SZ[k][j] = merk;
  }    
  if ( M->array.N )
  {
    for(k=0; k < M->rows; k++)
    {
      merk = M->array.N[k][i];
      M->array.N[k][i] = M->array.N[k][j];
      M->array.N[k][j] = merk;
    }
  }
  M->flags.Symmetric = FALSE;
}

/*}}}  */
/*{{{  row_add*/


/**************************************************************************\
@---------------------------------------------------------------------------
@ void row_add(M, i, j, fac)
@ matrix_TYP *M;
@ int i, j, fac;
@
@  adds fac times the i-th row the the j-th
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
row_add (matrix_TYP *M, int i, int j, int fac)
{
int k;

  for(k=0; k<M->cols; k++)
    M->array.SZ[j][k] += M->array.SZ[i][k] * fac;
  if ( M->prime != 0 )  {
    for(k=0; k<M->cols; k++) {
      M->array.SZ[j][k] %= M->prime;
      if ( M->array.SZ[j][k] < 0 ) M->array.SZ[j][k] = -M->array.SZ[j][k];
    }
  } else {
    if ( M->array.N != NULL )
      for(k=0; k<M->cols; k++)
        M->array.N[j][k] += M->array.N[i][k] * fac;
  }
}

/*}}}  */
/*{{{  col_add*/

/**************************************************************************\
@---------------------------------------------------------------------------
@ void row_col(M, i, j, fac)
@ matrix_TYP *M;
@ int i, j, fac;
@
@  adds fac times the i-th col the the j-th
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
col_add (matrix_TYP *M, int i, int j, int fac)
{
int k;

  for(k=0; k<M->rows; k++)
    M->array.SZ[k][j] += M->array.SZ[k][i] * fac;
  if ( M->prime != 0 )
  {
    for(k=0; k<M->rows; k++)
    {
      M->array.SZ[k][j] %= M->prime;
      if ( M->array.SZ[k][j] < 0 ) M->array.SZ[k][j] = -M->array.SZ[k][j];
    }
  }
  else
  {
    if ( M->array.N != NULL )
      for(k=0; k<M->rows; k++)
        M->array.N[k][j] += M->array.N[k][i] * fac;
  }
}

/*}}}  */
/*{{{  row_mul*/


/**************************************************************************\
@---------------------------------------------------------------------------
@ void row_mul(M, i, fac)
@ matrix_TYP *M;
@ int i, fac;
@
@  multiplies the i-th row with fac
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
row_mul (matrix_TYP *M, int i, int fac)
{
int k;

  for(k=0; k < M->cols; k++)
     M->array.SZ[i][k] *= fac;
  if ( M->prime != 0 )
  { 
    for(k=0; k < M->cols; k++)
    {
      M->array.SZ[i][k] %= M->prime;
      if ( M->array.SZ[i][k] < 0 ) M->array.SZ[i][k] = -M->array.SZ[i][k];
    }
  } 
  else
  {
    if ( M->array.N != NULL )
      for(k=0; k < M->cols; k++)
        M->array.N[i][k] *= fac;
  }
}

/*}}}  */
/*{{{  col_mul*/

/**************************************************************************\
@---------------------------------------------------------------------------
@ void col_mul(M, i, fac)
@ matrix_TYP *M;
@ int i, fac;
@
@  multiplies the i-th col with fac
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
col_mul (matrix_TYP *M, int i, int fac)
{
int k;

  for(k=0; k < M->rows; k++)
     M->array.SZ[k][i] *= fac;
  if (M->prime != 0 )
  {
    for(k=0; k < M->rows; k++)
    {
      M->array.SZ[k][i] %= M->prime;
      if ( M->array.SZ[k][i] < 0 ) M->array.SZ[k][i] = -M->array.SZ[k][i];
    }
  }
  else
  {
    if ( M->array.N != NULL )
      for(k=0; k < M->rows; k++)
        M->array.N[k][i] *= fac;
  }
}

/*}}}  */
