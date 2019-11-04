#include "typedef.h"
#include "gmp.h"
/* #include "gmp-impl.h" */
#include "longtools.h"
#include "matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  long_hnf.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ The function to the same as the functions in long_gauss.c,
@ they additionaly clear the columns upwards and transform 'Mat'
@ to a matrix in Hermite-normal-form
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ int long_row_hnf(Mat)
@ matrix_TYP *Mat;
@
@---------------------------------------------------------------------------
\**************************************************************************/
int 
long_row_hnf (matrix_TYP *Mat)
{
  int rang;
  MP_INT **M;
 
  M = matrix_to_MP_mat(Mat);
  rang = MP_hnf(M, Mat->rows, Mat->cols);
  write_MP_mat_to_matrix(Mat, M);
  free_MP_mat(M, Mat->rows, Mat->cols);
  free(M);
  return(rang);
}

int 
long_col_hnf (matrix_TYP *Mat)
{
  int i,
      j,
      rang;

  matrix_TYP *Mattr;

  Mattr = tr_pose(Mat);

  rang = long_row_hnf(Mattr);

  /* transpose the array */
  for (i=0;i<Mat->cols;i++)
     for (j=0;j<Mat->rows;j++)
        Mat->array.SZ[j][i] = Mattr->array.SZ[i][j];

  Check_mat(Mat);
  free_mat(Mattr);

  return rang;
}

/**************************************************************************\
@---------------------------------------------------------------------------
@ int long_row_trf_hnf(Mat, T)
@ matrix_TYP *Mat, *T;
@
@---------------------------------------------------------------------------
\**************************************************************************/
int 
long_row_trf_hnf (matrix_TYP *M, matrix_TYP *T)
{
  int rang;
  MP_INT **N, **S;

  N = matrix_to_MP_mat(M);
  S = init_MP_mat(M->rows, M->rows);
  rang = MP_trf_hnf(N, S, M->rows, M->cols);
  write_MP_mat_to_matrix(M, N);
  write_MP_mat_to_matrix(T, S);
  free_MP_mat(N, M->rows, M->cols); free_MP_mat(S, M->rows, M->rows);
  free(N); free(S);
  return(rang);
}




/**************************************************************************\
@---------------------------------------------------------------------------
@ int long_row_hnf_simultaneous(A, B)
@ matrix_TYP *A, *B;
@
@---------------------------------------------------------------------------
\**************************************************************************/
int 
long_row_hnf_simultaneous (matrix_TYP *A, matrix_TYP *B)
{
  int rang;
  MP_INT **MA, **MB;

  MA = matrix_to_MP_mat(A);
  MB = matrix_to_MP_mat(B);
  rang = MP_hnf_simultaneous(MA, A->rows, A->cols, MB, B->cols);
  write_MP_mat_to_matrix(A, MA);
  write_MP_mat_to_matrix(B, MB);
  free_MP_mat(MA, A->rows, A->cols); free_MP_mat(MB, B->rows, B->cols);
  free(MA); free(MB);
  return(rang);
}
