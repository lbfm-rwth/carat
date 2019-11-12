#include "typedef.h"
#include "utils.h"
#include "gmp.h"
#include "longtools.h"
#include "matrix.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: pair_red_inv.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *pair_red_inv(A, T)
@ matrix_TYP *A, *T;
@
@  calculates a matrices B and T such that  T A^{-1} T^{tr} = B 
@  is pair_reduced.
@  B is returned with B->kgv = 1 and the transformation is written to T.
@  The function used GNU MP to avoid overflow
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
pair_red_inv (matrix_TYP *A, matrix_TYP *T)
{
   MP_INT ***E, **MA, **MB, **Trf;
   MP_INT Ekgv;
   int Ecols;
   matrix_TYP *erg;
   int i,j;


   if(A->cols != A->rows)
   {
     printf("error: can't invert non-square matrix\n");
     exit(3);
   }
   MA = matrix_to_MP_mat(A);
   MB = (MP_INT **) xmalloc(A->cols *sizeof(MP_INT *));
   for(i=0;i<A->cols;i++)
   {
     MB[i] = (MP_INT *) xmalloc(A->cols *sizeof(MP_INT));
     for(j=0;j<A->cols;j++)
       mpz_init(&MB[i][j]);
     mpz_set_si(&MB[i][i], 1);
   }
   mpz_init(&Ekgv);
   E = MP_solve_mat(MA, A->rows, A->cols, MB, A->cols, &Ecols, &Ekgv);
   for(i=0;i<A->cols;i++)
   {
      for(j=0;j<A->cols;j++)
      { mpz_clear(&MA[i][j]); mpz_clear(&MB[i][j]);}
      free(MA[i]);
      free(MB[i]);
   }
   free(MA);
   free(MB);

   if(E[1] != NULL)
   {
     printf("Error: matrix in 'long_mat_inv' is singular\n");
     exit(3);
   }
   if(E[0] == NULL)
   {
     printf("Error in 'long_mat_inv'\n");
     exit(3);
   }
   Trf = (MP_INT **)xmalloc(A->cols *sizeof(MP_INT *));
   for(i=0;i<A->cols;i++)
   {
     Trf[i] = (MP_INT *)xmalloc(A->cols *sizeof(MP_INT));
     for(j=0;j<A->cols;j++)
        mpz_init(&Trf[i][j]);
     mpz_set_si(&Trf[i][i], 1);
   }
   MP_pair_red(E[0], Trf, A->cols);
   erg = MP_mat_to_matrix(E[0], A->cols, A->cols);
   write_MP_mat_to_matrix(T, Trf);
   for(i=0;i<A->cols;i++)
   {
      for(j=0;j<A->cols;j++)
      {
       mpz_clear(&E[0][i][j]);
       mpz_clear(&Trf[i][j]);
      }
      free(E[0][i]);
      free(Trf[i]);
   }
   free(E[0]);
   free(E);
   free(Trf);
   mpz_clear(&Ekgv);
   Check_mat(erg);
   return(erg);
}
