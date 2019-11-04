#include "typedef.h"
#include "gmp.h"
/* #include "gmp-impl.h" */
#include "longtools.h"
#include "matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  long_kernel_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/************************************************************************\
@ matrix_TYP *long_kernel_mat(A)
@ matrix_TYP *A;
@
@ long_kernel_mat(A) calculates Matrix X with AX = 0
@ the cols of X are a Z-basis of the solution space.
@ The functions uses GNU MP
\************************************************************************/
matrix_TYP *
long_kernel_mat (matrix_TYP *A)
{
   MP_INT ***E, **MA;
   MP_INT Ekgv;
   int Ecols, i,j;
   matrix_TYP *erg, *ergt;


   MA = matrix_to_MP_mat(A);
   mpz_init(&Ekgv);
   E = MP_solve_mat(MA, A->rows, A->cols, NULL, 0, &Ecols, &Ekgv);
   for(i=0;i<A->rows;i++)
   {
      for(j=0;j<A->cols;j++)
        mpz_clear(&MA[i][j]);
      free(MA[i]);
   }
   free(MA);


   if(E[1] != NULL)
   {
     erg = MP_mat_to_matrix(E[1], Ecols, A->cols);
     for(i=0;i<Ecols;i++)
     {
        for(j=0;j<A->cols;j++)
          mpz_clear(&E[1][i][j]);
        free(E[1][i]);
     }
     free(E[1]);
   }
   else
    erg = NULL;
   free(E);
   mpz_clear(&Ekgv);

   if (erg){
      ergt = tr_pose(erg);
   }
   else{
      ergt = NULL;
   }

   /* inserted the if 28/1/97 tilman, problems with a matrix having
      0 rows */
   if (erg != NULL){
      free_mat(erg);
   }
   return(ergt);
}
