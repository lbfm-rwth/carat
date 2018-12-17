#include "typedef.h"
#include "gmp.h"
#include "longtools.h"
/* #include "gmp-impl.h" */
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: MP_pair_red.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/



static void MP_two_reduce(G, T, i, j, n)
MP_INT **G, **T;
int i, j, n;
{
   int k;
   MP_INT f, f1, f2;
   MP_INT zwei;
  
   mpz_init(&f);
   mpz_init(&f1);
   mpz_init(&f2);
   mpz_init_set_si(&zwei, 2);
   mpz_div(&f1, &G[i][i], &zwei);
   if(mpz_cmp_si(&G[i][j], 0) > 0)
   {
     mpz_add(&f2, &G[i][j], &f1);
   }
   else
   {
     mpz_sub(&f2, &G[i][j], &f1);
   }
   mpz_div(&f, &f2, &G[i][i]);
   
   for(k=0;k<n;k++)
   {
       mpz_mul(&f1, &f, &G[i][k]);
       mpz_sub(&G[j][k], &G[j][k], &f1);
   }
   mpz_mul(&f1, &f, &G[j][i]);
   mpz_sub(&G[j][j], &G[j][j], &f1);
   for(k=0;k<n;k++)
   {
     if(k != j)
     {
       mpz_set(&G[k][j], &G[j][k]);
     }
   }
   for(k=0;k<n;k++)
   {
       mpz_mul(&f1, &f, &T[i][k]);
       mpz_sub(&T[j][k], &T[j][k], &f1);
   }
   mpz_clear(&f);
   mpz_clear(&f1);
   mpz_clear(&f2);
   mpz_clear(&zwei);
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ void MP_pair_red(G, T, n)
@ MP_INT **G, **T;
@ int   n;
@
@  Applies a pair reduction to the positive definite n x n - matrix G
@  and does the row operations simultaneous on T, i.e.
@  G is changed to T * G * T^{tr}
@  T must be initialized before calling the function.
@  It is possible to call the function with T = NULL.
@---------------------------------------------------------------------------
@
\**************************************************************************/
void MP_pair_red(G, T, n)
MP_INT **G, **T;
int   n;
{
   int i,j,k, p1, p2;
   int reduced, red2;
   MP_INT merk, s, a, zwei;

   mpz_init(&merk);
   mpz_init(&s);
   mpz_init(&a);
   mpz_init_set_si(&zwei, 2);
   reduced = FALSE;
   MP_reduction_sort(G,T,n);
   while(reduced == FALSE)
   {
     reduced = TRUE;
     for(i=0;i<n;i++)
      for(j=0;j<i;j++)
      {
         red2 = FALSE;
         while(red2 == FALSE && (mpz_cmp_si(&G[i][i], 0) != 0 || mpz_cmp_si(&G[j][j], 0) != 0))
         {
           if(mpz_cmp_si(&G[i][i], 0) > 0 && mpz_cmp(&G[i][i], &G[j][j])<= 0)
           {
             mpz_mul(&a, &zwei, &G[i][j]);
             mpz_abs(&merk, &a);
             if(mpz_cmp(&merk, &G[i][i]) > 0)
             {
               MP_two_reduce(G, T, i, j, n);
               reduced = FALSE;
             }
             else
              red2 = TRUE;
           }
           if(mpz_cmp_si(&G[j][j], 0) > 0 && mpz_cmp(&G[j][j], &G[i][i])<=0)
           {
             mpz_mul(&a, &zwei, &G[i][j]);
             mpz_abs(&merk, &a);
             if(mpz_cmp(&merk, &G[j][j]) > 0)
             {
               MP_two_reduce(G, T, j, i, n);
               reduced = FALSE;
             }
             else
              red2 = TRUE;
           }
           if(mpz_cmp_si(&G[i][i], 0) < 0 || mpz_cmp_si(&G[j][j], 0) < 0)
           {
             mpz_clear(&merk);
             mpz_clear(&s);
             mpz_clear(&a);
             mpz_clear(&zwei);
             return;
           }
         }
      }
      MP_reduction_sort(G, T, n);
      for(i=0;i<n;i++)
        for(j=i+1;j<n;j++)
        {
          if(mpz_cmp(&G[i][i], &G[j][j])<=0)
            { p1 = i; p2 = j;}
          else
            { p1 = j; p2 = i;}
          mpz_abs(&merk, &G[p1][p2]);
          if(mpz_cmp_si(&merk, 0) != 0)
            mpz_div(&s, &G[p1][p2], &merk);
          else
            mpz_set_si(&s, 0);
          mpz_mul(&merk, &merk, &zwei);
          if(mpz_cmp(&merk, &G[p1][p1]) == 0)
          {
            for(k=0; k<n;k++)
            {
              if(k != p2)
              {
                mpz_mul(&merk, &s, &G[p1][k]);
                mpz_sub(&G[p2][k], &G[p2][k], &merk);
              }
              mpz_mul(&merk, &s, &T[p1][k]);
              mpz_sub(&T[p2][k], &T[p2][k], &merk);
            }
            for(k=0;k<n;k++)
              mpz_set(&G[k][p2], &G[p2][k]);
            for(k=0;k<n;k++)
            {
              if(k != p1 && k != p2)
              {
                red2 = FALSE;
                while(red2 == FALSE && mpz_cmp_si(&G[k][k], 0) > 0 && mpz_cmp_si(&G[p2][p2], 0) > 0)
                {
                  mpz_abs(&a, &G[k][p2]);
                  mpz_mul(&a, &a, &zwei);
                  if(mpz_cmp(&a, &G[k][k])>0 || mpz_cmp(&a, &G[p2][p2])>0)
                  {
                     if(mpz_cmp(&a, &G[k][k])>0 && mpz_cmp(&G[k][k], &G[p2][p2])<0)
                      MP_two_reduce(G, T, k, p2, n);
                     else
                      MP_two_reduce(G, T, p2, k, n);
                     reduced = FALSE;
                  }
                  else
                   red2 = TRUE;
                }
                if(mpz_cmp_si(&G[k][k], 0) < 0 || mpz_cmp_si(&G[p2][p2], 0) < 0)
                {
                  mpz_clear(&merk);
                  mpz_clear(&s);
                  mpz_clear(&a);
                  mpz_clear(&zwei);
                  return;
                }
              }
            }
          }
          if(mpz_cmp_si(&G[i][i], 0) < 0 || mpz_cmp_si(&G[j][j], 0) < 0)
          {
            mpz_clear(&merk);
            mpz_clear(&s);
            mpz_clear(&a);
            mpz_clear(&zwei);
            return;
          }
        }
   }

   /* inserted 06/05/97 tilman */
   mpz_clear(&s);
   mpz_clear(&a);
   mpz_clear(&merk);
   mpz_clear(&zwei);

   return;

}
