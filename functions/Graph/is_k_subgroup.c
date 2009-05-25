/* last change: 13.02.2001 by Oliver Heidbuechel */


#include <typedef.h>
#include <presentation.h>
#include <matrix.h>
#include <bravais.h>
#include <base.h>
#include <datei.h>
#include <graph.h>
#include <gmp.h>
#include <zass.h>
#include <longtools.h>



/* -------------------------------------------------------------------------- */
/* Test if S is a proper k-subgroup of R                                      */
/* T(R) = Z^n                                                                 */
/* P = P(R) has to be P(S) (generators don't have to be equal)                */
/* (generating set without Id's from translations of the space group)         */
/* n: the last n matrices generate the translationlattice of S                */
/*    (we need this, if the presentation includes no information for this     */
/*     matrices)                                                              */
/* pres: presentation of P(R)                                                 */
/* -------------------------------------------------------------------------- */
boolean is_k_subgroup(bravais_TYP *S,
                      bravais_TYP *R,
                      bravais_TYP *P,
                      int n,
                      matrix_TYP *pres)
{
   matrix_TYP *coz_s_neu, *coz_r, *diff, **INV, *ID, **Id, *T;

   rational eins, minuseins;

   int i, j, k, d = S->gen[0]->rows - 1;

   word *relator;

   boolean FLAG = TRUE;

   bravais_TYP *S_neu;




   /* sublattice check for the last n generators */
   for (i = 1; i <= n; i++){
      Check_mat(S->gen[S->gen_no - i]);
      if (S->gen[S->gen_no - i]->kgv != 1)
         return(FALSE);
   }

   /* compare cocycles */
   coz_s_neu = sg(S, P);
   S_neu = extract_r(P, coz_s_neu);
   eins.z = eins.n = minuseins.n = 1; minuseins.z = -1;
   coz_r = extract_c(R);
   diff = mat_add(coz_s_neu, coz_r, eins, minuseins);
   Check_mat(diff);
   if (diff->kgv != 1){
      free_mat(coz_r);
      free_mat(coz_s_neu);
      free_mat(diff);
      free_bravais(S_neu);
      return(FALSE);
   }
   free_mat(diff);
   free_mat(coz_r);
   free_mat(coz_s_neu);

   /* compare translation lattices */
   INV = (matrix_TYP **)calloc(S_neu->gen_no, sizeof(matrix_TYP *));
   Id = (matrix_TYP **)calloc(pres->rows, sizeof(matrix_TYP *));
   relator = (word *)calloc(pres->rows, sizeof(word));
   for (i = 0; i < pres->rows; i++){
      matrix_2_word(pres, relator + i, i);
   }
   for (i = 0; i < pres->rows; i++){
      Id[i] = matrizen_in_word(S_neu->gen, INV, relator[i]);
   }
   T = init_mat(d, n + pres->rows, "");
   for (i = 0; i < n; i++){
      for (j = 0; j < d; j++)
         T->array.SZ[j][i] = S->gen[S->gen_no - 1 - i]->array.SZ[j][d];
   }
   for (i = 0; i < pres->rows && FLAG; i++){
      if (Id[i]->kgv != 1)
         FLAG = FALSE;
      for (j = 0; j < d && FLAG; j++){    /* rows */
         for (k = 0; k < d && FLAG; k++){ /* colums */
            if (j == k){
               if (Id[i]->array.SZ[j][k] != 1)
                  FLAG = FALSE;
            }
            else{
               if (Id[i]->array.SZ[j][k] != 0)
                  FLAG = FALSE;
            }
         }
      }
      for (j = 0; j < d && FLAG; j++){ /* last column */
         T->array.SZ[j][i + n] = Id[i]->array.SZ[j][d];
      }
   }
   ID = init_mat(d, d, "1");
   if (FLAG == TRUE){
      long_col_hnf(T);
      real_mat(T, T->rows, T->rows);

      if (cmp_mat(T, ID) != 0)
         FLAG = TRUE;
      else
         FLAG = FALSE;
   }

   /* clean */
   free_mat(ID);
   free_mat(T);
   for (i = 0; i < S_neu->gen_no; i++){
      if (INV[i] != NULL)
         free_mat(INV[i]);
   }
   free(INV);
   free_bravais(S_neu);
   for (i = 0; i < pres->rows; i++)
      free_mat(Id[i]);
   free(Id);
   for (i = 0; i < pres->rows; i++)
      wordfree(relator + i);
   free(relator);

   return(FLAG);
}
