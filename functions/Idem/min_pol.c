#include "typedef.h"
#include "matrix.h"
#include "longtools.h"


/***************************************************************************
@
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: min_pol.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
****************************************************************************/


/***************************************************************************
@
@---------------------------------------------------------------------------
@
@ matrix_TYP *min_pol(matrix_TYP *A)
@
@ Let A be an integral NON ZERO matrix. Then the function returns an
@ integral 1xn matrix B with the property:
@ \sum_i=0^n B->array.SZ[0][i] * A^i =0,
@ and n is minimal.
@
@ SIDEEFECTS: the matrix A is checked via Check_mat !!!
@---------------------------------------------------------------------------
@
****************************************************************************/
matrix_TYP *min_pol(matrix_TYP *A)
{
   int i,
       j,
       k,
       flag,
       rank,
       grad;     /* actualy the degree +1 */

   matrix_TYP **potenzen,
               *erg;        /* will hold the result */

   MP_INT  **lines,
           **trf,
             kgv;

   /* check some trivialities */
   Check_mat(A);
   if (A->flags.Integral == FALSE || A->cols != A->rows){
      fprintf(stderr,"error in min_pol, A should be integral and quadratic\n");
      exit(3);
   }

   /* prepare all variables */
   potenzen = (matrix_TYP **) calloc(1 + A->cols , sizeof(matrix_TYP *));
   potenzen[0] = init_mat(A->cols,A->cols,"1");
   potenzen[1] = copy_mat(A);
   lines = init_MP_mat(A->cols+1,A->cols*A->cols);
   trf = init_MP_mat(A->cols+1,A->cols+1);
   for (i=0;i<=A->cols;i++) mpz_set_ui(trf[i]+i,1);
   mpz_init(&kgv);

   /* set the first two lines */
   for (i=0;i<2;i++)
      for (j=0;j<A->cols;j++)
         for (k=0;k<A->cols;k++)
            mpz_set_si(lines[i]+j*A->cols+k,potenzen[i]->array.SZ[j][k]);

   /* the calculation starts */
   grad = 2;
   do{
      if (INFO_LEVEL & 4){
         dump_MP_mat(lines,grad,A->cols*A->cols,"lines");
      }
      rank = MP_row_gauss_simultaneous(lines,grad,A->cols*A->cols,trf,grad);

      /* avoid getting trouble by big numbers, ie. the powers aren't really
         powers because of a lack of precision */
      if (rank > A->cols){
         fprintf(stderr,"error in min_pol because of a lack of precision\n");
         exit(3);
      }

      if (rank == grad){
         potenzen[grad] = mat_mul(potenzen[grad-1],A);
         for (j=0;j<A->cols;j++)
            for (k=0;k<A->cols;k++)
               mpz_set_si(lines[grad]+j*A->cols+k,
                     potenzen[grad]->array.SZ[j][k]);
         grad++;
         flag = FALSE;
      }
      else{
         flag = TRUE;
      }
   } while(!flag);

   /* we got the wanted result in trf */
   erg = init_mat(1,grad,"");
   for (i=0;i<grad;i++) erg->array.SZ[0][i] = mpz_get_si(&trf[grad-1][i]);

   /* "normalize" polynomial, ie. get the last coeff positive */
   if (erg->array.SZ[0][grad-1] < 0)
      for (i=0;i<grad;i++) erg->array.SZ[0][i] *= -1 ;

   /* clean up memory */
   for (i=0;i<=A->cols;i++) if (potenzen[i] != NULL) free_mat(potenzen[i]);
   free_MP_mat(trf,A->cols+1,A->cols+1);
   free(trf);
   free(potenzen); 
   free_MP_mat(lines,A->cols+1,A->cols*A->cols);
   free(lines);
   mpz_clear(&kgv);

   return erg;
} /* min_pol(....) */
