#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <longtools.h>
#include <presentation.h>


matrix_TYP *mapped_word(int *w,
                        matrix_TYP **A,
                        matrix_TYP **AINV)
{

   int i;

   matrix_TYP *M;

   if (w[0] == 0)
      return init_mat(A[0]->cols,A[0]->cols,"i1");

   if (w[1] < 0) {
     if (AINV[-w[1]-1] == NULL)
       AINV[-w[1]-1] = long_mat_inv(A[-w[1]-1]);
     M = copy_mat(AINV[-w[1]-1]);
   }
   else{
     M = copy_mat(A[w[1]-1]);
   }

   for (i=2;i<=w[0];i++){
      if (w[i] < 0){
         if (AINV[-w[i]-1] == NULL)
            AINV[-w[i]-1] = long_mat_inv(A[-w[i]-1]);
         mat_muleq(M,AINV[-w[i]-1]);
      }
      else{
         mat_muleq(M,A[w[i]-1]);
      }
   }

   return M;

} /* mapped_word(...) */


