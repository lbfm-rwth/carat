#include "typedef.h"
#include "gmp.h"
/* #include "gmp-impl.h" */
#include "longtools.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: long_elt.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *long_elt_mat(Mat, left_trans)
@ matrix_TYP *Mat, *left_trans;
@
@ calculates the elementary divisors of the matrix Mat, t.m.
@ calulates a 'diagonal' matrix D such that
@     L * Mat * R = D for some matrices L and R
@ 'left_trans' has to be a matrix of size Mat->rows x Mat->rows
@ it is changed to L * left_trans.
@ It is possible to start the function with left_trans = NULL.
@
@ The calculations are done using multiple precision integers.
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *long_elt_mat(Mat, left_trans)
matrix_TYP *Mat, *left_trans;
{

   int i,j,k, step, stepclear;
   int cols, rows, rang;
   MP_INT **trf, **M, **Mt, *merkpointer;
   matrix_TYP *erg;
   int t_option;
   MP_INT a1, a2, x1, x2, y1, y2, merk, g, f;
 
   cols = Mat->cols;
   rows = Mat->rows;
   if(left_trans != NULL)
   {
     t_option = TRUE;
     if(left_trans->cols != rows)
     {
       printf("error in 'long_elt_mat':\n");
       printf("the matrix 'left_trans' has to be a %dx%d-matrix, but has %d  columns\n", rows, cols, left_trans->cols);
       exit(3);
     }
     if(left_trans->cols != left_trans->rows)
     {
       printf("error: matrix 'left_trans' in 'long_elt_met' has to be a square matrix\n");
       exit(3);
     }
   }
   else
    t_option = FALSE;
   mpz_init(&a1), mpz_init(&a2);
   mpz_init(&x1), mpz_init(&x2);
   mpz_init(&y1), mpz_init(&y2);
   mpz_init(&f), mpz_init(&g);
   mpz_init(&merk);
   /***************************************************************\
   | Set Mt= Mat^{tr} transform Mt to Hermite normal form.
   \***************************************************************/
   if((Mt = (MP_INT **)malloc(cols *sizeof(MP_INT *))) == NULL)
   {
     printf("malloc of 'Mt' in 'long_elt_mat' failed\n");
     exit(2);
   }
   for(i=0;i<cols;i++)
   {
     if((Mt[i] = (MP_INT *)malloc(rows *sizeof(MP_INT))) == NULL)
     {
       printf("malloc of 'Mt[%d]' in 'long_elt_mat' failed\n", i);
       exit(2);
     }
     for(j=0;j<rows;j++)
       mpz_init_set_si(&Mt[i][j], Mat->array.SZ[j][i]);
   }
   rang = MP_hnf(Mt, cols, rows);
   /***************************************************************\
   | Set trf =  left_tans
   \***************************************************************/
   if(t_option == TRUE)
   {
     if((trf = (MP_INT **)malloc(rows *sizeof(MP_INT *))) == NULL)
     {
       printf("malloc of 'trf' in 'long_elt_mat' failed\n");
       exit(2);
     }
     for(i=0;i<rows;i++)
     {
       if((trf[i] = (MP_INT *)malloc(rows *sizeof(MP_INT))) == NULL)
       {
         printf("malloc of 'trf[%d]' in 'long_elt_mat' failed\n", i);
         exit(2);
       }
       for(j=0;j<rows;j++)
         mpz_init_set_si(&trf[i][j], left_trans->array.SZ[i][j]);
     }
   }
   /***************************************************************\
   | Set M= Mt^{tr} transform Mt to Hermite normal form.
   \***************************************************************/
   if((M = (MP_INT **)malloc(rows *sizeof(MP_INT *))) == NULL)
   {
     printf("malloc of 'M' in 'long_elt_mat' failed\n");
     exit(2);
   }
   for(i=0;i<rows;i++)
   {
     if((M[i] = (MP_INT *)malloc(cols *sizeof(MP_INT))) == NULL)
     {
       printf("malloc of 'M[%d]' in 'long_elt_mat' failed\n", i);
       exit(2);
     }
     for(j=0;j<cols;j++)
       mpz_init_set(&M[i][j], &Mt[j][i]);
   }
   if(t_option == TRUE)
     rang = MP_hnf_simultaneous(M, rows, cols, trf, rows);
   else
     rang = MP_hnf(M, rows, cols);
   /***************************************************************\
   | Clear the space allocated for 'Mt'
   \***************************************************************/
   for(i=0;i<cols;i++)
   {
     for(j=0;j<rows;j++)
       mpz_clear(&Mt[i][j]);
     free(Mt[i]);
   }
   free(Mt);
   /***************************************************************\
   | Now the elementary divisor algorithm starts
   \***************************************************************/
   for(step = 0; step < rang;step++)
   {

      do
      {
         /* output for debugging
         dump_MP_mat(M,rang,rang,"M"); */

         /*------------------------------------------------------*\
         | Clear the 'step'^th row of M
         \*------------------------------------------------------*/
          for(i=step;i<rang && mpz_cmp_si(&M[step][i], 0) == 0; i++);
          if(i!=step)
          {
            for(j=step;j<rang;j++)
            {
              mpz_set(&merk, &M[j][step]);
              mpz_set(&M[j][step], &M[j][i]);
              mpz_set(&M[j][i], &merk);
            }
          }
          if(mpz_cmp_si(&M[step][step], 0) < 0)
          {
            for(i=step;i<rows;i++)
              mpz_neg(&M[i][step], &M[i][step]);
          }
          for(i=step+1;i<rang;i++)
          {
            if(mpz_cmp_si(&M[step][i], 0) != 0)
            {
              if(mpz_cmp_si(&M[step][step], 1) == 0)
              {
                mpz_set(&f, &M[step][i]);
                for(j=step+1;j<rang;j++)
                {
                  mpz_mul(&merk, &M[j][step], &f);
                  mpz_sub(&M[j][i], &M[j][i], &merk);
                }
                mpz_set_si(&M[step][i], 0);

                /* commented the next 8 lines out: tilman 5/12/96
                if(t_option == TRUE)
                {
                   for(j=0;j<rows;j++)
                   {
                     mpz_mul(&merk, &trf[j][step], &f);
                     mpz_sub(&trf[j][i], &trf[j][i], &merk);
                   }
                }*/

              }
              else
              {
                 mpz_gcdext(&g, &a1, &a2, &M[step][step], &M[step][i]);
                 mpz_div(&x1, &M[step][i], &g);
                 mpz_div(&x2, &M[step][step], &g);
                 mpz_neg(&x2, &x2);

                 for(j=step+1; j<rang;j++)
                 {
                    mpz_mul(&y1, &a1, &M[j][step]);
                    mpz_mul(&merk, &a2, &M[j][i]);
                    mpz_add(&y1, &y1, &merk);

                    mpz_mul(&y2, &x1, &M[j][step]);
                    mpz_mul(&merk, &x2, &M[j][i]);
                    mpz_add(&y2, &y2, &merk);

                    mpz_set(&M[j][step], &y1);
                    mpz_set(&M[j][i], &y2);
                 }
                 mpz_set(&M[step][step], &g);
                 mpz_set_si(&M[step][i], 0);

                 /* commented the next 16 lines out: tilman 5/12/96
                 if(t_option == TRUE)
                 {
                   for(j=0; j<rows;j++)
                   {
                      mpz_mul(&y1, &a1, &trf[j][step]);
                      mpz_mul(&merk, &a2, &trf[j][i]);
                      mpz_add(&y1, &y1, &merk);

                      mpz_mul(&y2, &x1, &trf[j][step]);
                      mpz_mul(&merk, &x2, &trf[j][i]);
                      mpz_add(&y2, &y2, &merk);

                      mpz_set(&trf[j][step], &y1);
                      mpz_set(&trf[j][i], &y2);
                   }
                 }*/

              }
            }
          }
         /*------------------------------------------------------*\
         | Clear the 'step'^th column of M
         \*------------------------------------------------------*/
          for(i=step;i<rang && mpz_cmp_si(&M[i][step], 0) == 0; i++);
          if(i!=step)
          {
            merkpointer = M[step];
            M[step] = M[i];
            M[i] = merkpointer;
            if(t_option == TRUE)
            {
              merkpointer = trf[step];
              trf[step] = trf[i];
              trf[i] = merkpointer;
            }
          }
          if(mpz_cmp_si(&M[step][step], 0) < 0)
          {
            for(i=step;i<rows;i++)
              mpz_neg(&M[i][step], &M[i][step]);
          }
          for(i=step+1;i<rang;i++)
          {
            if(mpz_cmp_si(&M[i][step], 0) != 0)
            {
              if(mpz_cmp_si(&M[step][step], 1) == 0)
              {
                mpz_set(&f, &M[i][step]);
                for(j=step+1;j<rang;j++)
                {
                  mpz_mul(&merk, &M[step][j], &f);
                  mpz_sub(&M[i][j], &M[i][j], &merk);
                }
                mpz_set_si(&M[i][step], 0);

                /* inserted the next 7 lines: tilman 5/12/96 */
                if (t_option){
                   for(j=0;j<rows;j++)
                   {
                     mpz_mul(&merk, &trf[step][j], &f);
                     mpz_sub(&trf[i][j], &trf[i][j], &merk);
                   }
                }

              }
              else
              {
                 /* changed on 8/1/97 tilman from
                 mpz_gcdext(&g, &a1, &a2, &M[step][step], &M[i][step]);
                 to : */
                 mpz_gcdext(&g, &a1, &a2, &M[step][step], &M[step][i]);

                 mpz_div(&x1, &M[i][step], &g);
                 mpz_div(&x2, &M[step][step], &g);
                 mpz_neg(&x2, &x2);

                 for(j=step+1; j<rang;j++)
                 {
                    mpz_mul(&y1, &a1, &M[step][j]);
                    mpz_mul(&merk, &a2, &M[i][j]);
                    mpz_add(&y1, &y1, &merk);

                    mpz_mul(&y2, &x1, &M[step][j]);
                    mpz_mul(&merk, &x2, &M[i][j]);
                    mpz_add(&y2, &y2, &merk);

                    mpz_set(&M[step][j], &y1);
                    mpz_set(&M[i][j], &y2);
                 }

                 /* inserted the next 15 lines: tilman 5/12/96 */
                 if (t_option){
                    for(j=0; j<rows;j++)
                    {
                       mpz_mul(&y1, &a1, &trf[step][j]);
                       mpz_mul(&merk, &a2, &trf[i][j]);
                       mpz_add(&y1, &y1, &merk);

                       mpz_mul(&y2, &x1, &trf[step][j]);
                       mpz_mul(&merk, &x2, &trf[i][j]);
                       mpz_add(&y2, &y2, &merk);

                       mpz_set(&trf[step][j], &y1);
                       mpz_set(&trf[i][j], &y2);
                    }
                 }

                 mpz_set(&M[step][step], &g);
                 mpz_set_si(&M[i][step], 0);
              }
            }
          }
         /*------------------------------------------------------*\
         | Check whether the column and row are both clear
         \*------------------------------------------------------*/
         stepclear = TRUE;
         for(i=step+1;i<cols && stepclear == TRUE;i++)
         {
           if(mpz_cmp_si(&M[step][i], 0) != 0)
             stepclear = FALSE;
         }
         for(i=step+1;i<rows && stepclear == TRUE;i++)
         {
           if(mpz_cmp_si(&M[i][step], 0) != 0)
             stepclear = FALSE;
         }
      }while(stepclear == FALSE);
      if(mpz_cmp_si(&M[step][step], 0) < 0)
         mpz_neg(&M[step][step], &M[step][step]);
   }
   /*******************************************************************\
   | Now M is diagonal.
   | Now one has transform M such that M[i][i] divides M[i+1][i+1]
   \*******************************************************************/
   for(step = 0;step < rang; step ++)
   {
      /*******************************************************************\
      | Search for the minimal entry of M[step][step],...,M[rang][rang]
      | and swap it to M[step][step]
      \*******************************************************************/
      j = step;
      for(i=step+1; i<rang;i++)
      {
         if(mpz_cmp(&M[i][i], &M[j][j]) < 0)
            j = i;
      }
      if(j != step)
      {
        mpz_set(&merk, &M[step][step]);
        mpz_set(&M[step][step], &M[j][j]);
        mpz_set(&M[j][j], &merk);
        if(t_option != 0)
        {
           merkpointer = trf[step];
           trf[step] = trf[j];
           trf[j] = merkpointer;
        }
      }
      for(i=step+1;i<rang;i++)
      {
         mpz_mod(&merk, &M[i][i], &M[step][step]);
         if(mpz_cmp_si(&merk, 0) != 0)
         {
           if(t_option == TRUE)
           {
            mpz_gcdext(&g, &a1, &a2, &M[step][step], &M[i][i]);
            for(j=0;j<rows;j++)
               mpz_add(&trf[step][j], &trf[step][j], &trf[i][j]);

            mpz_div(&f, &M[i][i], &g);
            mpz_mul(&f, &f, &a2);
            for(j=0;j<rows;j++)
            {
               mpz_mul(&merk, &f, &trf[step][j]);
               mpz_sub(&trf[i][j], &trf[i][j], &merk);
            }
            mpz_div(&merk, &M[i][i], &g);
            mpz_mul(&M[i][i], &merk, &M[step][step]);

            /* changed 28/2/97 tilman from
            mpz_div(&M[step][step], &M[step][step], &g);
            to: */
            mpz_set(&M[step][step],&g);
           }
           else
           {
             mpz_gcd(&g, &M[step][step], &M[i][i]);
             mpz_div(&merk, &M[i][i], &g);
             mpz_mul(&M[i][i], &merk, &M[step][step]);

             /* changed 28/2/97 tilman from
             mpz_div(&M[step][step], &M[step][step], &g);
             to: */
             mpz_set(&M[step][step],&g);
           }
         }
      }
   }
   erg = MP_mat_to_matrix(M, rows, cols);
   if(t_option == TRUE)
   {
     for(i=0;i<rows;i++)
       for(j=0;j<rows;j++)
       {
         if(abs(trf[i][j]._mp_size) > 1)
         {
           printf("Error: Integer overflow in 'long_mat_elt'\n");
           exit(3);
         }
         left_trans->array.SZ[i][j] = mpz_get_si(&trf[i][j]);
       }
   }
   for(i=0;i<rows;i++)
   {
      if(t_option == TRUE)
      {
        for(j=0;j<rows;j++)
          mpz_clear(&trf[i][j]);
        free(trf[i]);
      }
      for(j=0;j<cols;j++)
        mpz_clear(&M[i][j]);
      free(M[i]);
   }
   free(M);
   if(t_option == TRUE)
     free(trf);
   merkpointer = NULL;
   mpz_clear(&a1), mpz_clear(&a2);
   mpz_clear(&x1), mpz_clear(&x2);
   mpz_clear(&y1), mpz_clear(&y2);
   mpz_clear(&f), mpz_clear(&g);
   mpz_clear(&merk);
   Check_mat(erg);
   return(erg);
}
