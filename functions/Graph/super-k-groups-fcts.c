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
#include <tsubgroups.h>



extern int OANZ;
extern matrix_TYP **OMAT;
extern int OFLAG;



/* -------------------------------------------------------------------- */
/* R: spacegroup in standard affine form without translations           */
/* calculate all minimal k-supergroups of R                             */
/* T(R) has to be Z^n                                                   */
/* P is the point group of R with correct normalizer                    */
/* return the number of these groups via anz                            */
/* if bahnflag representatives under the operation of N_Aff_n (R)       */
/* are returned                                                         */
/* debugflag == TRUE: do paranoia test                                  */
/* -------------------------------------------------------------------- */
bravais_TYP **min_k_super(bravais_TYP *P,
                          bravais_TYP *R,
			  matrix_TYP *pres,
                          int *anz,
			  boolean bahnflag,
                          boolean debugflag,
			  int **length)
{
   bravais_TYP **S,
                *TR;

   matrix_TYP *F, *std, *tmp2,
              **gitter, **allgitter,
	      **N;

   int p, i, j, Nanz, anzalt,
       *list = 0, *smallest = 0, *kopiert;


   anz[0] = 0;
   OFLAG = TRUE;
   OANZ = 0;
   OMAT = (matrix_TYP **)malloc(1000 * sizeof(matrix_TYP *));
   F = init_mat(P->dim, P->dim, "1");

   /* calculate the superlattices */
   for (p = 2; p < 100 ; p++){
      if (P->divisors[p]){
         TR = tr_bravais(P, 1, 0);
         memset(TR->divisors, 0, 100 * sizeof(int));
         TR->divisors[p] = 1;

	 /* TR ist changed in ZZ */
	 ZZ(TR, F, TR->divisors, NULL, "rbgl1", 0, 0, 0);
	 free_bravais(TR);

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
   free_mat(F);

   /* calculate the correct lattices */
   allgitter = (matrix_TYP **)calloc(anz[0], sizeof(matrix_TYP *));
   for (i = 0; i < anz[0]; i++){
      allgitter[i] = mat_inv(OMAT[i]);
      long_col_hnf(allgitter[i]);
      free_mat(OMAT[i]);
   }
   free(OMAT);

   if (bahnflag){
      /* calculate the point group of the affine normalizer of R */
      cen_to_norm(P);
      N = PoaN(R, P, pres, &Nanz);

      /* calculate orbit on lattices */
      list = (int *)calloc(anz[0], sizeof(int));
      length[0] = (int *)calloc(anz[0] + 1, sizeof(int));
      smallest = (int *)calloc(anz[0] + 1, sizeof(int));
      anzalt = anz[0];
      anz[0] = orbit_on_lattices(allgitter, anzalt, N, Nanz, list,
                                 length[0], smallest, NULL);

      /* get representatives */
      gitter = (matrix_TYP **)calloc(anz[0], sizeof(matrix_TYP *));
      kopiert = (int *)calloc(anzalt, sizeof(int));
      for (i = 1; i <= anz[0]; i++){
         gitter[i - 1] = allgitter[smallest[i]];
	 kopiert[smallest[i]] = 1;
      }

      /* clean */
      for (i = 0; i < Nanz; i++){
         free_mat(N[i]);
      }
      for (i = 0; i < anzalt; i++){
         if (kopiert[i] == 0)
	    free_mat(allgitter[i]);
      }
      free(kopiert);
      free(N);
      free(allgitter);
   }
   else{
      gitter = allgitter;
   }

   /* calculate the groups */
   S = (bravais_TYP **)calloc(anz[0], sizeof(bravais_TYP *));
   for (i = 0; i < anz[0]; i++){
      S[i] = copy_bravais(R);
      plus_translationen(S[i], gitter[i]);

      /* paranoia test */
      if (debugflag){
         std = copy_mat(gitter[i]);
         long_col_hnf(std);
         for (j = 0; j < P->gen_no; j++){
            tmp2 = mat_mul(P->gen[j], gitter[i]);
            long_col_hnf(tmp2);
            if (cmp_mat(std, tmp2)){
               fprintf(stderr, "Not G-invariant!!!\n");
               exit(2);
            }
            free_mat(tmp2);
         }
         free_mat(std);
      }
   }

   /* clean */
   for (i = 0; i < anz[0]; i++){
      free_mat(gitter[i]);
   }
   free(gitter);
   if (bahnflag){
      free(list);
      free(smallest);
   }

   return(S);
}
























