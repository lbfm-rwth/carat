#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <longtools.h>
#include <presentation.h>



/* -------------------------------------------------------------------------- */
/* tests if M == Id                                                           */
/* -------------------------------------------------------------------------- */
static boolean equal_id(matrix_TYP *M)
{
   int i, j;

   if (M->rows != M->cols)
      return(FALSE);

   for (i = 0; i < M->rows; i++){
      for (j = 0; j < M->cols; j++){
         if (i == j){
            if (M->array.SZ[i][j] != 1)
               return(FALSE);
         }
         else{
            if (M->array.SZ[i][j] != 0)
               return(FALSE);
         }
      }
   }
   return(TRUE);
}



/* -------------------------------------------------------------------------- */
/* transform a matrix A to standard form                                      */
/* -------------------------------------------------------------------------- */
void standard_form(matrix_TYP *A,
                   matrix_TYP *D,
                   int first)
{
   int i, j;

   for (i = 0; i < A->rows; i++){
      for (j = 0; j < A->cols; j++){
         A->array.SZ[j][i] %= D->array.SZ[j + first][j + first];
         if (A->array.SZ[j][i] < 0)
            A->array.SZ[j][i] += D->array.SZ[j + first][j + first];
      }
   }
}



/* -------------------------------------------------------------------------- */
/* invert A : C_a x C_b x ... -> C_a x C_b x ...                              */
/* -------------------------------------------------------------------------- */
matrix_TYP *graph_mat_inv(matrix_TYP *A,
                          matrix_TYP *D,
                          int first)
{
   matrix_TYP *B, *C;

   C = copy_mat(A);
   B = init_mat(A->rows, A->cols, "1");

   while (equal_id(C) != TRUE){
      free_mat(B);
      B = C;
      C = NULL;
      C = mat_mul(B, A);
      standard_form(C, D, first);
   }
   free_mat(C);
   return(B);
}



/* -------------------------------------------------------------------------- */
/* mapped word for A[i] : C_a x C_b x ... -> C_a x C_b x ...                  */
/* -------------------------------------------------------------------------- */
matrix_TYP *graph_mapped_word(int *w,
                              matrix_TYP **A,
                              matrix_TYP **AINV,
                              matrix_TYP *D)
{

   int i, first;

   matrix_TYP *M;


   for (first = 0; first < D->cols && D->array.SZ[first][first] == 1; first++);

   if (w[0] == 0)
      return init_mat(A[0]->cols,A[0]->cols,"i1");

   if (w[1] < 0) {
     if (AINV[-w[1]-1] == NULL)
       AINV[-w[1]-1] = graph_mat_inv(A[-w[1]-1], D, first);
     M = copy_mat(AINV[-w[1]-1]);
   }
   else{
     M = copy_mat(A[w[1]-1]);
   }

   for (i=2;i<=w[0];i++){
      if (w[i] < 0){
         if (AINV[-w[i]-1] == NULL)
            AINV[-w[i]-1] = graph_mat_inv(A[-w[i]-1], D, first);
         mat_muleq(M,AINV[-w[i]-1]);
      }
      else{
         mat_muleq(M,A[w[i]-1]);
      }
      standard_form(M, D, first);
   }

   return M;

} /* graph_mapped_word(...) */


