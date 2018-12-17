/* last change: 07.02.2001 by Oliver Heidbuechel */


#include <typedef.h>
#include <presentation.h>
#include <matrix.h>
#include <bravais.h>
#include <base.h>
#include <datei.h>
#include <graph.h>
#include <gmp.h>
#include <zass.h>
#include <tools.h>
#include <longtools.h>
#include <voronoi.h>
#include "../functions/ZZ/ZZ_P.h"
#include <ZZ.h>


extern int OANZ;
extern matrix_TYP **OMAT;
extern int OFLAG;
/* boolean GRAPH = FALSE; */




/* -------------------------------------------------------------------- */
/* P < GL_n(Z) finite                                                   */
/* Returns the P-invariant maximal sublattices of Z^n                   */
/* (matrices in long_col_hnf - Form)				        */
/* with p_i-power-index, where p_i is a prime given in P->divisors.     */
/* Returns the number of these sublattices via anz.                     */
/* set trivialflag[x] = TRUE iff lattices[x] = p_i * Z^n                */
/* -------------------------------------------------------------------- */
matrix_TYP **max_sublattices(bravais_TYP *P,
                             int *anz,
                             int **trivialflag,
                             boolean debugflag)
{
   matrix_TYP *F, *tmp, *std, **lattices;

   int p, X[100], i, j;

   bravais_TYP *PP;


   anz[0] = 0;
   OFLAG = TRUE;
   OANZ = 0;
   memcpy(X, P->divisors, 100 * sizeof(int));
   OMAT = (matrix_TYP **)calloc(1000, sizeof(matrix_TYP *));
   trivialflag[0] = (boolean *)calloc(1000, sizeof(boolean));
   F = init_mat(P->dim, P->dim, "1");

   for (p = 2; p < 100 ; p++){
      if (X[p]){
         PP = copy_bravais(P);
         memset(PP->divisors, 0, 100 * sizeof(int));
         PP->divisors[p] = 1;

	 /* PP is changed in ZZ */
	 ZZ(PP, F, PP->divisors, NULL, "rbgl1", 0, 0, 0);
         free_bravais(PP);

         if (OANZ == anz[0]){
            /* trivial lattice isn't returned */
            OMAT[OANZ] = init_mat(P->dim, P->dim, "");
            for(i = 0; i < OMAT[OANZ]->rows; i++){
               OMAT[OANZ]->array.SZ[i][i] = p;
            }
            Check_mat(OMAT[OANZ]);
            trivialflag[0][OANZ] = TRUE;
            OANZ++;
         }
         else{
            for (i = anz[0]; i < OANZ; i++){
               ZZ_transpose_array(OMAT[i]->array.SZ, OMAT[i]->cols);
               long_col_hnf(OMAT[i]);
               trivialflag[0][i] = FALSE;
            }
         }
         cleanup_prime();
      }
      anz[0] = OANZ;
   }

   OANZ = 0;
   OFLAG = FALSE;
   free_mat(F);
   lattices = (matrix_TYP **)calloc(anz[0], sizeof(matrix_TYP *));

   for (i = 0; i < anz[0]; i++){

      /* paranoia test */
      if (debugflag){
         std = copy_mat(OMAT[i]);
         long_col_hnf(std);
         for (j = 0; j < P->gen_no; j++){
            tmp = mat_mul(P->gen[j], OMAT[i]);
            long_col_hnf(tmp);
            if (cmp_mat(std, tmp)){
               fprintf(stderr, "Not G-invariant!!!\n");
               exit(2);
            }
            free_mat(tmp);
         }
         free_mat(std);
      }

      lattices[i] = copy_mat(OMAT[i]);
      free_mat(OMAT[i]);
   }
   free(OMAT);
   memcpy(P->divisors, X, 100 * sizeof(int));

   return(lattices);
}

