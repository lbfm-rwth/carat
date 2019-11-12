#include "typedef.h"
#include "utils.h"
#include "gmp.h"
/* #include "gmp-impl.h" */
#include "longtools.h"
#include "matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: long_solve_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/************************************************************************\
@-------------------------------------------------------------------------
@ matrix_TYP **long_solve_mat(A, B)
@ matrix_TYP *A, *B;
@
@ long_solve_mat(A,B) calculates Matrix X[0] with AX[0] = B, and
@ a matrix X[1] with MX[1] = 0, such that
@ the cols of X[1] are a Z-basis of the solution space.
@-------------------------------------------------------------------------
\************************************************************************/
matrix_TYP **
long_solve_mat (matrix_TYP *A, matrix_TYP *B)
{
   MP_INT ***E, **MA, **MB;
   MP_INT Ekgv;
   int Ecols, i,j;
   matrix_TYP **erg;
   matrix_TYP *et;
   int Bcols;

   MA = matrix_to_MP_mat(A);
   if(B != NULL)
   {
     MB = matrix_to_MP_mat(B);
     Bcols = B->cols;
   }
   else
   {
     MB = NULL; Bcols = 0;
   }
   mpz_init(&Ekgv);
   E = MP_solve_mat(MA, A->rows, A->cols, MB, Bcols, &Ecols, &Ekgv);
   for(i=0;i<A->rows;i++)
   {
      for(j=0;j<A->cols;j++)
        mpz_clear(&MA[i][j]);
      free(MA[i]);
   }
   free(MA);
   if(B != NULL)
   {
     /* changed tilman on 21/2/97 from:
     for(i=0;i<A->rows;i++)
     to: */
     for(i=0;i<B->rows;i++)
     {
        for(j=0;j<Bcols;j++)
          mpz_clear(&MB[i][j]);
        free(MB[i]);
     }
     free(MB);
   }


   erg = (matrix_TYP **)xmalloc(2 *sizeof(matrix_TYP *));
   if(E[0] != NULL)
   {
     erg[0] = MP_mat_to_matrix(E[0], A->cols, Bcols);
     erg[0]->flags.Integral = FALSE;
     if(abs(Ekgv._mp_size) > 1)
     {
       printf("Error: Integer overflow in 'MP_mat_to_matrix'\n");
       exit(3);
     }
     erg[0]->kgv = mpz_get_si(&Ekgv);
     if(A->kgv != 1)
     {
       for(i=0;i<A->cols;i++)
         for(j=0;j<Bcols;j++)
           erg[0]->array.SZ[i][j] *= A->kgv;
     }
     if(B != NULL && B->kgv != 1)
       erg[0]->kgv *= B->kgv;
     Check_mat(erg[0]);
     for(i=0;i<A->cols;i++)
     {
        for(j=0;j<Bcols;j++)
          mpz_clear(&E[0][i][j]);
        free(E[0][i]);
     }
     free(E[0]);
   }
   else
     erg[0] = NULL;
   if(E[1] != NULL)
   {
     et = MP_mat_to_matrix(E[1], Ecols, A->cols);
     for(i=0;i<Ecols;i++)
     {
        for(j=0;j<A->cols;j++)
          mpz_clear(&E[1][i][j]);
        free(E[1][i]);
     }
     free(E[1]);
     erg[1] = tr_pose(et);
     free_mat(et);
   }
   else
    erg[1] = NULL;
   free(E);
   mpz_clear(&Ekgv);
   return(erg);
}
