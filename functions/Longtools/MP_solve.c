#include "typedef.h"
#include "gmp.h"
#include "longtools.h"
/* #include "gmp-impl.h" */

extern int INFO_LEVEL;

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  MP_solve.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/************************************************************************\
@-------------------------------------------------------------------------
@ MP_INT ***MP_solve_mat(M, rows, cols, B, Bcols, X1cols, X0kgv)
@ MP_INT **M, **B, *X0kgv;
@ int rows, cols, Bcols, *X1cols;
@
@ MP_solve_mat(M) calculates Matrix X[0] with MX[0] = B, and
@ a matrix X[1] with MX[1] = 0, such that
@ the cols of X[1] are a Z-basis of the solution space.
@ The number of the cols of X[1] are returned via the pointer 'X1cols'.
@ CAUTION: the function changes the values of M;
@-------------------------------------------------------------------------
@
\************************************************************************/
MP_INT ***MP_solve_mat(M, rows, cols, B, Bcols, X1cols, X0kgv)
MP_INT **M, **B, *X0kgv;
int rows, cols, Bcols, *X1cols;
{
   MP_INT **Mt, **Trf, ***erg, *Y, *Z;
   MP_INT kgv, altkgv, merk, merk1, g;
   MP_INT zaehler, nenner, b;
   int i,j,k,l,rang, n;
   int tester, dim;
   int waste;

   extern matrix_TYP *copy_mat();
   extern matrix_TYP *tr_pose();
   extern matrix_TYP *init_mat();
   extern matrix_TYP *mat_mul();

   mpz_init(&zaehler); mpz_init(&nenner); mpz_init(&b);
   if(B == NULL || Bcols == 0)
     rang = MP_row_gauss(M, rows, cols);
   else
     rang = MP_row_gauss_simultaneous(M, rows, cols, B, Bcols);
   if((Mt = (MP_INT **)malloc(cols *sizeof(MP_INT *))) == 0)
   {
     printf("malloc of 'Mt' in 'MP_solve_mat' failed\n");
     exit(2);
   }
   for(i=0;i<cols;i++)
   {
     if((Mt[i] = (MP_INT *)malloc(rows *sizeof(MP_INT))) == 0)
     {
       printf("malloc of 'Mt[%d]' in 'MP_solve_mat' failed\n", i);
       exit(2);
     }
   }
   for(i=0;i<rows;i++)
     for(j=0;j<cols;j++)
       mpz_init_set(&Mt[j][i], &M[i][j]);
   if((Trf = (MP_INT **)malloc(cols *sizeof(MP_INT *))) == 0)
   {
     printf("malloc of 'Trf' in 'MP_solve_mat' failed\n");
     exit(2);
   }
   for(i=0;i<cols;i++)
   {
     if((Trf[i] = (MP_INT *)malloc(cols *sizeof(MP_INT))) == 0)
     {
       printf("malloc of 'Trf[%d]' in 'MP_solve_mat' failed\n", i);
       exit(2);
     }
     for(j=0;j<cols;j++)
        mpz_init(&Trf[i][j]);
     mpz_set_si(&Trf[i][i], 1);
   }
   mpz_init_set_si(&kgv,1);

   /* transform the transposed Mt of M to be gauss reduced */
   rang = MP_trf_gauss(Mt, Trf, cols, rows);

   /* now we are able to predict the dimension of the homogenous solution */
   dim = cols - rang;
   *X1cols = dim;

   if((erg = (MP_INT ***)malloc(2 *sizeof(MP_INT **))) == 0)
   {
     printf("malloc of 'erg' in 'MP_solve_mat' failed\n");
     exit(2);
   }

   {
      /*********************************************************************\
      | calculate a solution of the homogenous equation
      | Write this solution to erg[1]
      \*********************************************************************/
     if(dim == 0)
     {
        erg[1] = NULL;
     }
     else
     {
       if((erg[1] = (MP_INT  **)malloc(dim *sizeof(MP_INT *))) == 0)
       {
         printf("malloc of 'erg[1]' in 'MP_solve_mat' failed\n");
         exit(2);
       }

       for(i=0;i<dim;i++)
       {
         if((erg[1][i] = (MP_INT  *)malloc(cols *sizeof(MP_INT))) == 0)
         {
           printf("malloc of 'erg[1][%d]' in 'MP_solve_mat' failed\n", i);
           exit(2);
         }
       }
       for(i=0;i<dim;i++)
          for(j=0;j<cols;j++)
             mpz_init_set(&erg[1][i][j], &Trf[i+rang][j]);
     }
     if(dim > 1)
       MP_hnf(erg[1], dim, cols);
   }
   if(B != NULL && Bcols != 0)
   {

     /* test whether there is an solution at all */
     /* does work this way because M and B were gauss reduced
        simultaneusly */
     tester = TRUE;
     for(i=rang;i<rows && tester == TRUE;i++)
       for(j=0;j<Bcols && tester == TRUE;j++)
       {
          if(mpz_cmp_si(&B[i][j], 0) != 0)
            tester = FALSE;
       }

     if(tester == FALSE)
       erg[0] = NULL;
     else
     {
      /*********************************************************************\
      | calculate a solution of the inhomogenous equation
      | Write this solution to erg[0]
      \*********************************************************************/
      if ((erg[0] = (MP_INT **)malloc(cols *sizeof(MP_INT *))) == NULL)
      {
         printf("malloc of 'erg[0]' in 'MP_solve_mat' failed\n");
         exit(2);
       }

       for(i=0;i<cols;i++)
       {
         if((erg[0][i] = (MP_INT  *)malloc(Bcols *sizeof(MP_INT))) == 0)
         {
           printf("malloc of 'erg[0][%d]' in 'MP_solve_mat' failed\n", i);
           exit(2);
         }
         for(j=0;j<Bcols;j++)
           mpz_init(&erg[0][i][j]);
       }
      if ((Y = (MP_INT *)malloc(cols *sizeof(MP_INT))) == NULL)
      {
        printf("malloc of 'Y' in 'MP_solve_mat' failed\n");
        exit(2);
      }
      for(i=0;i<cols;i++)
        mpz_init(&Y[i]);
      if ((Z = (MP_INT *)malloc(cols *sizeof(MP_INT))) == NULL)
      {
        printf("malloc of 'Z' in 'MP_solve_mat' failed\n");
        exit(2);
      }
      for(i=0;i<cols;i++)
        mpz_init(&Z[i]);
      mpz_init_set_si(&altkgv,1);
      mpz_init(&merk); mpz_init(&merk1);
      mpz_init(&g);
      for(i=0;i<Bcols;i++)
      {

        /* output for debugging */
        if (INFO_LEVEL & 4){
           dump_MP_mat(M,rows,cols,"M");
           dump_MP_mat(Mt,cols,rows,"Mt");
           dump_MP_mat(B,rows,Bcols,"B");
        }

        /* changed 30/1/97 tilman from
        for(j=0;j<cols && j< rows;j++)
        to : */
        for(j=0;j<cols && j< rang;j++)
        {

          if (INFO_LEVEL & 4){
             dump_MP_mat(&Y,1,cols,"Y");
             dump_MP_mat(&Z,1,cols,"Z");
          }

          mpz_set_si(&Y[j], 0);
          mpz_set_si(&Z[j], 1);
          for(k=0;k<j;k++)
          {
             mpz_mul(&merk1, &Y[k], &Mt[k][j]);
             mpz_mul(&merk1, &merk1, &Z[j]);
             mpz_mul(&Y[j], &Y[j], &Z[k]);
             mpz_add(&Y[j], &Y[j], &merk1);
             mpz_mul(&Z[j], &Z[j], &Z[k]);

             mpz_gcd(&g, &Y[j], &Z[j]);
             mpz_div(&Y[j], &Y[j], &g);
             mpz_div(&Z[j], &Z[j], &g);
          }
          mpz_mul(&merk, &B[j][i], &Z[j]);
          mpz_sub(&Y[j], &merk, &Y[j]);
          mpz_mul(&Z[j], &Z[j], &Mt[j][j]);

          mpz_gcd(&g, &Y[j], &Z[j]);

          mpz_div(&Y[j], &Y[j], &g);
          mpz_div(&Z[j], &Z[j], &g);

          /* normalising the sign of Z[j] */
          if(mpz_cmp_si(&Z[j], 0) < 0){
            mpz_neg(&Y[j], &Y[j]);
            mpz_neg(&Z[j], &Z[j]);
          }
        }

        /* changed 30/1/97 from:
        for(j=rows;j<cols;j++)
        to :*/
        for(j=rang;j<cols;j++)
        {
          mpz_set_si(&Y[j], 0);
          mpz_set_si(&Z[j], 1);
        }
        /* calculate the lcm of kgv and Z[j], 1<=j<=cols */
        mpz_set(&merk, &kgv);
        for(j=0;j<cols && j < rows;j++)
        {
             mpz_gcd(&g, &Z[j], &merk);
             mpz_div(&merk1, &Z[j], &g);
             mpz_mul(&merk, &merk, &merk1);
        }
        for(j=0;j<cols;j++)
        {
          mpz_div(&g, &merk, &Z[j]);
          mpz_mul(&Y[j], &Y[j], &g);
        }
        if(mpz_cmp(&merk, &kgv) != 0)
        {
          mpz_div(&g, &merk, &kgv);
          mpz_set(&kgv, &merk);
          for(j=0;j<Bcols;j++)
          {
            if(j != i)
            {
              for(k=0;k<cols;k++)
                mpz_mul(&erg[0][k][j], &erg[0][k][j], &g);
            }
          }
        }
        for(j=0;j<cols;j++)
        {
           mpz_set_si(&erg[0][j][i], 0);
           for(k=0;k<cols;k++)
           {
              mpz_mul(&merk, &Y[k], &Trf[k][j]);
              mpz_add(&erg[0][j][i], &erg[0][j][i], &merk);
           }
        }
      }

      /* output for debugging */
      if (INFO_LEVEL & 4){
         fprintf(stderr,"kgv = ");
         mpz_out_str(stderr,10,&kgv);
         fprintf(stderr,"\n");
         dump_MP_mat(erg[0],cols,Bcols,"erg[0]");
      }

      mpz_clear(&altkgv);
      mpz_clear(&merk); mpz_clear(&merk1);
      mpz_clear(&g);
      for(i=0;i<cols;i++)
        mpz_clear(&Y[i]);
      free(Y);
      for(i=0;i<cols;i++)
        mpz_clear(&Z[i]);
      free(Z);
   }
   mpz_set(X0kgv, &kgv);
  }
   else
     erg[0] = NULL;

   for(i=0;i<cols;i++)
   {
     for(j=0;j<cols;j++)
       mpz_clear(&Trf[i][j]);
     free(Trf[i]);
   }
   free(Trf);
   for(i=0;i<rows;i++)
     for(j=0;j<cols;j++)
       mpz_clear(&Mt[j][i]);
   for(i=0;i<cols;i++)
      free(Mt[i]);
   free(Mt);

   mpz_clear(&zaehler);
   mpz_clear(&nenner);
   mpz_clear(&b);
   mpz_clear(&kgv);

   return(erg);
}
