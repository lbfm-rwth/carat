#include "typedef.h"
#include "gmp.h"
#include "getput.h"
#include "zass.h"
#include "matrix.h"

void convert_cocycle_to_column(matrix_TYP **Y,
                               int number,
                               int dim,
                               int gen_no)
{
   matrix_TYP *A;

   int i,
       j,
       k,
       kgv;


   for (i=0;i<number;i++){

       if (dim != Y[i]->rows || gen_no != Y[i]->cols){
          fprintf(stderr,"error in convert_cocycle_to_column\n");
          exit(3);
       }
       Check_mat(Y[i]);

       A = init_mat(dim*gen_no,1,"i");

       for (j=0;j<Y[i]->rows;j++){
          for (k=0;k<Y[i]->cols;k++){
             A->array.SZ[j+k*dim][0] = Y[i]->array.SZ[j][k];
          }
       }

       kgv = Y[i]->kgv;
       free_mat(Y[i]);
       Y[i] = A;
       Y[i]->kgv  = kgv;
       Check_mat(Y[i]);

   }

   return;

}

