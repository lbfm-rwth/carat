#include "typedef.h"
#include "matrix.h"
#include "longtools.h"
#include "getput.h"
#include "tools.h"
#include "zass.h"


matrix_TYP **cong_solve(matrix_TYP *A)
/* loest das system linearer congruenzen A x = 0 ( mod Z^n) fuer
   x in Q^n.
   Die loesungen stehen in cong_solve[0], die zugehoerigen nenner
   in cong_solve[1], die transformationsmatrix in cong_solve[2] */
{
    matrix_TYP **loesungen,
               *umform,
               *A_tr;
    long i,
         j,
         found,
         hilf;

    if (INFO_LEVEL & 4){
      printf("entering cong_solve\n");
    }

    loesungen = (matrix_TYP **) malloc(4 * sizeof(matrix_TYP*));

    umform = scalar(A->cols,1);

    A_tr = tr_pose(A);
    loesungen[1] = long_elt_mat(umform,A_tr,NULL);
    loesungen[1]->kgv = A->kgv;
    loesungen[1]->flags.Integral = A->flags.Integral;
    Check_mat(loesungen[1]);
    loesungen[2] = tr_pose(umform);
    loesungen[3] = copy_mat(loesungen[1]);
    real_mat(loesungen[3], loesungen[3]->rows, loesungen[3]->rows);

    if (INFO_LEVEL & 4){
       put_mat(loesungen[1],NULL,"loesungen[1]",2);
       put_mat(umform,NULL,"umform",2); 
    }

    real_mat(loesungen[1],A->cols,A->cols);

    if (INFO_LEVEL & 4){
       put_mat(A,NULL,"A",2);
       put_mat(loesungen[1],NULL,"loesungen[1]",2);
       put_mat(umform,NULL,"umform",2);
    }

    loesungen[0] = init_mat(A->cols,1,"k");
    found = 0;

    for (i=0;i<loesungen[1]->cols;i++){
      if (((loesungen[1]->array.SZ[i][i])!=1)  &&
          ((loesungen[1]->array.SZ[i][i])!=(-1))){

         /* die matrix der loesungen wird um eine spalte groesser */
         found++;
         real_mat(loesungen[0],loesungen[0]->rows,found);

         /* speichern der neu gefundenen loesung */
         hilf = loesungen[0]->cols-1;
         for (j=0;j<loesungen[0]->rows;j++){
            loesungen[0]->array.SZ[j][hilf] =
                          umform->array.SZ[i][j];
         }

         /* speichern des dazugehoerigen nenners */
         loesungen[1]->array.SZ[found-1][found-1] = loesungen[1]->array.SZ[i][i];
         /* loesungen[1]->array.SZ[hilf][hilf] = loesungen[1]->array.SZ[i][i]; */
      } 
    }

    free_mat(umform);
    free_mat(A_tr);

    if (found == 0)
       real_mat(loesungen[0],loesungen[0]->rows,found);
    real_mat(loesungen[1],loesungen[1]->rows,loesungen[0]->cols);
    real_mat(loesungen[1],loesungen[0]->cols,loesungen[0]->cols);

    if (INFO_LEVEL & 4){
       put_mat(loesungen[0],NULL,"loesungen[0]",2);
       put_mat(loesungen[1],NULL,"loesungen[1]",2);
    }

    Check_mat(loesungen[0]);
    Check_mat(loesungen[1]);

    return loesungen;
}  /* cong_solve */
