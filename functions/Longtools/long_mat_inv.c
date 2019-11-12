#include "typedef.h"
#include "utils.h"
#include "gmp.h"
/* #include "gmp-impl.h"*/
#include "longtools.h"
#include "matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: long_mat_inv.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/************************************************************************\
@ matrix_TYP *long_mat_inv(A)
@ matrix_TYP *A;
@
@ long_mat_inv(A) calculates the matrix A^{-1} using GNU MP
\************************************************************************/
matrix_TYP *
long_mat_inv (matrix_TYP *A)
{
   MP_INT ***E, **MA, **MB;
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
   erg = MP_mat_to_matrix(E[0], A->cols, A->cols);
   if(abs(Ekgv._mp_size) > 1)
   {
     printf("Error: Integer overflow in 'MP_mat_to_matrix'\n");
     exit(3);
   }
   erg->kgv = mpz_get_si(&Ekgv);
   for(i=0;i<A->cols;i++)
   {
      for(j=0;j<A->cols;j++)
       mpz_clear(&E[0][i][j]);
      free(E[0][i]);
   }
   free(E[0]);
   free(E);
   mpz_clear(&Ekgv);
   if(A->kgv != 1)
   {
     for(i=0;i<A->cols;i++)
       for(j=0;j<A->cols;j++)
         erg->array.SZ[i][j] *= A->kgv;
   }
   Check_mat(erg);
   return(erg);
}
