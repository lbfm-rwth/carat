#include "typedef.h"
#include "gmp.h"
/* #include "gmp-impl.h" */
#include "longtools.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: long_gauss.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ int long_row_gauss(Mat)
@ matrix_TYP *Mat;
@
@  applies an integral gaussian algorithm to the rows of mat.
@  the rank is returned
@  the function uses GNU MP to avoid overflow
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
long_row_gauss (matrix_TYP *Mat)
{
  int rang;
  MP_INT **M;
 
  M = matrix_to_MP_mat(Mat);
  rang = MP_row_gauss(M, Mat->rows, Mat->cols);
  write_MP_mat_to_matrix(Mat, M);
  free_MP_mat(M, Mat->rows, Mat->cols);
  free(M);
  return(rang);
}

/************************************************************************
@
@------------------------------------------------------------------------
@
@ int long_row_basis(Mat,int flag)
@
@ returns the rank of the INTEGRAL matrix Mat.
@ The matrix Mat is changed in such a way that the first rank rows
@ of it form an integral basis of the Lattices spanned by the rows
@ of the matrix Mat.
@
@------------------------------------------------------------------------
@
*************************************************************************/
int 
long_row_basis (matrix_TYP *Mat, int flag)
{
  int rang,
      cols = Mat->cols,
      i,
      j,
      k;

  MP_INT   merk,
         **M,
         **M3,
         **M4,
         **T1;

 
  M = matrix_to_MP_mat(Mat);
  rang = MP_row_gauss(M, Mat->rows, Mat->cols);

  /* write this matrix back and check if it fits into int's */
  k = FALSE;
  for (i=0;i<Mat->rows;i++)
     for (j=0;j<Mat->cols;j++)
       if (abs(M[i][j]._mp_size)>1)
          k = TRUE;
       else
          Mat->array.SZ[i][j] = mpz_get_si(&M[i][j]);

  if (k || flag){
     /***************************************************************\
     | make a pair_reduction to the rows of M if needed
     \***************************************************************/
     mpz_init_set_ui(&merk,0);
     M3 = init_MP_mat(rang, rang);
     for(i=0;i<rang;i++)
       for(j=0;j<=i;j++)
       {
         for(k=0;k<cols;k++)
         {
          mpz_mul(&merk, &M[i][k], &M[j][k]);
          mpz_add(&M3[i][j], &merk, &M3[i][j]);
         }
       }
     for(i=0;i<rang;i++)
       for(j=0;j<i;j++)
         mpz_set(&M3[j][i], &M3[i][j]);

     /* output for debugging purposes
      dump_MP_mat(M3,rows,rows,"M3"); */

      T1 = init_MP_mat(rang, rang);
      for(i=0;i<rang;i++)
        mpz_set_si(&T1[i][i], 1);
      MP_pair_red(M3, T1, rang);

      /* output for debugging purposes
      dump_MP_mat(T1,rows,rows,"T1"); */

      M4 = init_MP_mat(Mat->rows, cols);
      for(i=0;i<rang;i++)
       for(j=0;j<cols;j++)
         for(k=0;k<rang;k++)
         {
            mpz_mul(&merk, &T1[i][k], &M[k][j]);
            mpz_add(&M4[i][j], &M4[i][j], &merk);
         }
   
      /* output for debugging purposes
      dump_MP_mat(M4,rows,rows,"M4"); */

     write_MP_mat_to_matrix(Mat, M4);

     free_MP_mat(M, Mat->rows, Mat->cols);
     free(M);
     free_MP_mat(M3,rang,rang);
     free(M3);
     free_MP_mat(M4,Mat->rows,cols);
     free(M4);
     free_MP_mat(T1, rang,rang);
     free(T1);
     mpz_clear(&merk);
  }
  else{
     free_MP_mat(M, Mat->rows, Mat->cols);
     free(M);
  }

  return(rang);
}

/**************************************************************************\
@---------------------------------------------------------------------------
@ int long_row_trf_gauss(Mat,T)
@ matrix_TYP *Mat, *T;
@
@  applies an integral gaussian algorithm to the rows of mat.
@  the rank is returned
@  the function uses GNU MP to avoid overflow.
@  The tranformation is written to the matrix T.
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
long_row_trf_gauss (matrix_TYP *M, matrix_TYP *T)
{
  int rang;
  MP_INT **N, **S;

  N = matrix_to_MP_mat(M);
  S = init_MP_mat(M->rows, M->rows);
  rang = MP_trf_gauss(N, S, M->rows, M->cols);
  write_MP_mat_to_matrix(M, N);
  write_MP_mat_to_matrix(T, S);
  free_MP_mat(N, M->rows, M->cols); free_MP_mat(S,M->rows, M->rows);
  free(N); free(S);
  return(rang);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ int long_row_gauss_simultaneous(A, B)
@ matrix_TYP *A, *B;
@
@  Applies an integral gaussian algorithm to the rows of mat.
@  The rank is returned
@  The function uses GNU MP to avoid overflow.
@  The tranformations are applied simultaneously to the matrix B.
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
long_row_gauss_simultaneous (matrix_TYP *A, matrix_TYP *B)
{
  int rang;
  MP_INT **MA, **MB;

  MA = matrix_to_MP_mat(A);
  MB = matrix_to_MP_mat(B);
  rang = MP_row_gauss_simultaneous(MA, A->rows, A->cols, MB, B->cols);
  write_MP_mat_to_matrix(A, MA);
  write_MP_mat_to_matrix(B, MB);
  free_MP_mat(MA, A->rows, A->cols); free_MP_mat(MB, B->rows, B->cols);
  free(MA); free(MB);
  return(rang);
}
