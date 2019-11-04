#include"typedef.h"
#include"bravais.h"
#include"matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: normlin.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *normlin(Fo, N, fdim)
@ matrix_TYP **Fo, *N;
@ int fdim;
@
@ normlin calculates a matrix fdim x fdim - matrix X such that
@ X[i], the i-th row of X, has the property
@ X[i][0] * Fo[0] + ... + X[i][fdim-1] * Fo[fdim-1] = N^{tr} * Fo[i] * N
@ That means X describes the matrix of the linear action of N on
@ the vectorspace generated by the Fo[i] by F->N^{tr}F N with respect to
@ the basis Fo[0],..,Fo[fdim-1].
@ CAUTION: The matrix describes the action on rows !
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
normlin (matrix_TYP **Fo, matrix_TYP *N, int fdim)
{
  int i,j,k, dim;
  matrix_TYP *A, *waste, *erg, *Ntr;

  dim = Fo[0]->cols;
  if(N->cols != dim || N->rows != dim)
  {
    printf("wrong dimension of 'N' in 'normlin'\n");
    exit(3);
  }
  erg = init_mat(fdim, fdim, "");
  Ntr = tr_pose(N);
  for(i=0;i<fdim;i++)
  {
    waste =  mat_mul(Ntr, Fo[i]);
    A = mat_mul(waste, N);
    form_to_vec(erg->array.SZ[i], A, Fo, fdim, &j);
    if(j == -1)
    {
      for(k=0;k<fdim;k++)
        erg->array.SZ[i][k] = -erg->array.SZ[i][k];
      j = 1;
    }
    if(j != 1){
     printf("Error in 'normlin': no Z-basis of formspace\n");
     exit(3);
    }
    free_mat(A);
    free_mat(waste);
  }
  free_mat(Ntr);
  return(erg);
}
