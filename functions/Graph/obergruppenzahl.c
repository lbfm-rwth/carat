#include <typedef.h>
#include <base.h>
#include <matrix.h>
#include <longtools.h>
#include <zass.h>
#include <tools.h>
#include <graph.h>





/* ----------------------------------------------------------------------------- */




int obergruppenzahl(matrix_TYP *L,
                    matrix_TYP **Norm,
                    matrix_TYP **StabStdCoz,
		    int Stab_anz,
		    int *wort)
{
   matrix_TYP *N, *NL, **orbit, *tmp;

   int i, j, k, counter;


   /* Normalisatorelement, welches zum Standardvertreter konjugiert */
   if (wort != NULL){
      N = init_mat(Norm[0]->rows, Norm[0]->rows, "1");
      for (i = 1; i <= wort[0]; i++){
         mat_muleq(N, Norm[wort[i]]);
      }
      NL = mat_mul(N, L);
      free_mat(N);
   }
   else{
      NL = copy_mat(L);
   }

   /* transformiere Gitter */
   long_col_hnf(NL);

   /* Berechne Bahn */
   counter = 1;
   orbit = (matrix_TYP **)calloc(32, sizeof(matrix_TYP *));
   orbit[0] = NL;
   for (i = 0; i < counter; i++){
      for (j = 0; j < Stab_anz; j++){
         tmp = mat_mul(StabStdCoz[j], orbit[i]);
	 long_col_hnf(tmp);
	 for (k = 0; k < counter; k++){
	    if (cmp_mat(tmp, orbit[k]) == 0)
	       break;
	 }
	 if (k == counter){
	    orbit[counter] = tmp;
	    counter++;
	    if (counter % 32 == 0)
	       orbit = (matrix_TYP **)realloc(orbit, (counter + 32) * sizeof(matrix_TYP *));
	 }
	 else{
	    free_mat(tmp);
	 }
      }
   }

   /* Aufraeumen */
   for (i = 0; i < counter; i++)
      free_mat(orbit[i]);
   free(orbit);

   return(counter);
}






