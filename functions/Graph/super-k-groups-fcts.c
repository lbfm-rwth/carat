/* last change: 19.01.01 by Oliver Heidbuechel */


#include <ZZ.h>
#include <typedef.h>
#include <presentation.h>
#include <matrix.h>
#include <bravais.h>
#include <base.h>
#include <datei.h>
#include <longtools.h>
#include <graph.h>


extern int OANZ;
extern matrix_TYP **OMAT;
extern int OFLAG;
boolean GRAPH = FALSE;


/* -------------------------------------------------------------------- */
/* calculate all minimal k-supergroups of R                             */
/* T(R) has to be Z^n                                                   */
/* P is the point group of R                                            */
/* return the number of these groups via anz                            */
/* -------------------------------------------------------------------- */
bravais_TYP **min_k_super(bravais_TYP *P,
                          bravais_TYP *R,
                          int *anz,
                          boolean debugflag)
{
   bravais_TYP **S,
                *TR;

   matrix_TYP *F, *tmp, *std, *tmp2;

   int p, X[100], i, j;


   anz[0] = 0;
   OFLAG = TRUE;
   OANZ = 0;
   memcpy(X, P->divisors, 100 * sizeof(int));
   OMAT = (matrix_TYP **)malloc(1000 * sizeof(matrix_TYP *));
   TR = tr_bravais(P, 1, 0);
   F = init_mat(P->dim, P->dim, "1");

   /* calculate the superlattices */
   for (p = 2; p < 100 ; p++){
      if (X[p]){
         memset(TR->divisors, 0, 100 * sizeof(int));
         TR->divisors[p] = 1;

         ZZ(TR, F, P->divisors, NULL, "rbgl1", 0, 0, 0);

         if (OANZ == anz[0]){
            OMAT[OANZ] = init_mat(P->dim, P->dim, "");
            for(i = 0; i < OMAT[OANZ]->rows; i++){
               OMAT[OANZ]->array.SZ[i][i] = p;
            }
            Check_mat(OMAT[OANZ]);
            OANZ++;
         }
         cleanup_prime();
      }
      anz[0] = OANZ;
   }

   OANZ = 0;
   OFLAG = FALSE;
   free_bravais(TR);
   free_mat(F);

   /* calculate the groups */
   S = (bravais_TYP **)calloc(anz[0], sizeof(bravais_TYP *));
   for (i = 0; i < anz[0]; i++){
      S[i] = copy_bravais(R);
      tmp = mat_inv(OMAT[i]);

      /* paranoia test */
      if (debugflag){
         std = copy_mat(tmp);
         long_col_hnf(std);
         for (j = 0; j < P->gen_no; j++){
            tmp2 = mat_mul(P->gen[j], tmp);
            long_col_hnf(tmp2);
            if (cmp_mat(std, tmp2)){
               fprintf(stderr, "Not G-invariant!!!\n");
               exit(2);
            }
            free_mat(tmp2);
         }
         free_mat(std);
      }

      plus_translationen(S[i], tmp);
      free_mat(OMAT[i]);
      free_mat(tmp);
   }
   free(OMAT);

   return(S);
}
























