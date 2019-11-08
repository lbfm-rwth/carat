#include "typedef.h"
#include "idem.h"
#include "matrix.h"
#include "bravais.h"
#include "longtools.h"
#include "tools.h"

matrix_TYP *almost_decomposable_lattice(bravais_TYP *G)
{

   int i,
       j,
       k,
       dimc,
       dimcc,
       col,
       den,
       IDEM_NO;

   matrix_TYP *id,
              *F,
             **IDEM_SPACES,
             **IDEM,
              *new_base,
              *tmp;


   id = init_mat(G->dim,G->dim,"1");

   /* we need a G-invariant, positive definite form for various reasons */
   F = rform(G->gen,G->gen_no,id,101);


   /* get the idempotents of the group */
   IDEM = idempotente(G->gen,G->gen_no,F,&IDEM_NO,&dimc,&dimcc);
   IDEM_SPACES = (matrix_TYP **) malloc(IDEM_NO * sizeof(matrix_TYP *));
   den = 1; for (i=0;i<IDEM_NO;i++) den = KGV(den,IDEM[i]->kgv);
   for (i=0;i<IDEM_NO;i++){
      /* tmp = tr_pose(IDEM[i]);
      IDEM_SPACES[i] = long_rein_mat(tmp);
      free_mat(tmp); */
      tmp = tr_pose(IDEM[i]);
      tmp->kgv = 1; iscal_mul(tmp,den / IDEM[i]->kgv);
      j = long_row_gauss(tmp);
      real_mat(tmp,j,tmp->cols);
      IDEM_SPACES[i] = tmp;
   }

   new_base = init_mat(G->dim,G->dim,"0");
   col = -1;
   for (i=0;i<IDEM_NO;i++){
      for (j=0;j<IDEM_SPACES[i]->rows;j++){
         col++;
         for (k=0;k<G->dim;k++){
            new_base->array.SZ[k][col] = IDEM_SPACES[i]->array.SZ[j][k];
         }
      }
   }

   Check_mat(new_base);
   tmp = mat_inv(new_base); tmp->kgv=1; Check_mat(tmp);
   free_mat(new_base);
   for (i=0;i<IDEM_NO;i++){
      free_mat(IDEM_SPACES[i]);
   }
   for (i=0;i<IDEM_NO+dimc+dimcc;i++){
      free_mat(IDEM[i]);
   }
   free(IDEM_SPACES);
   free(IDEM);
   free_mat(F);
   free_mat(id);

   return tmp;

}
