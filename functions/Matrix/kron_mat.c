#include"typedef.h"
#include"matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: kron_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *kron_mat(A,B)
@ matrix_TYP *A, *B;
@
@ calculates the Kronecker product of the matrices A and B.
@---------------------------------------------------------------------------
@
\**************************************************************************/

matrix_TYP *
kron_mat (matrix_TYP *A, matrix_TYP *B)
{

  int i,j,k,l,r,c;
  int rows, cols;
  matrix_TYP *C;

  rows = A->rows *B->rows;
  cols = A->cols *B->cols;
  C = init_mat(rows,cols, "");
  for(i=0; i<A->rows;i++)
    for(j=0;j<A->cols;j++)
    {
       /* changed 15/4/97 tilman
       r = i*A->rows; */
       r = i*B->rows;
       for(k=0;k<B->rows;k++)
       {
        /* changed 15/4/97 tilman
        c = j*A->cols; */
        c = j*B->cols;
        for(l=0;l<B->cols;l++)
        {
          C->array.SZ[r][c] = A->array.SZ[i][j] * B->array.SZ[k][l];
          c++;
        }
        r++;
       }
    }
  return(C);
}