#include "typedef.h"
#include "gmp.h"
/* #include "gmp-impl.h" */

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: MP_hnf.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
@
@ The algorithms do the same as the algorithms in the
@ file MP_gauss.c, the only difference the the matrix is in
@ Hermite-normal-form after the algorithm
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ int MP_trf_hnf(M, Trf, rows, cols)
@ MP_INT **M, **Trf;
@ int rows, cols;
@
@---------------------------------------------------------------------------
\**************************************************************************/
/************************************************************************\
| The Matrix M is transformed to a matrix M' such that
| M' = Trf * M is Gauss-reduced, with Trf integral with determinante +-1.
| CAUTION: The entries of the matrix M are changed.
\************************************************************************/

int 
MP_trf_hnf (MP_INT **M, MP_INT **Trf, int rows, int cols)
{
  int i,j;
  int n;
  int step;
  MP_INT a1,a2,gcd, *v, f;
  int tester;
  int spos = 0;
  MP_INT x1,x2,y1,y2;
  MP_INT Mi, Ms;
  int rang = 0;

  n = rows;
  mpz_init(&f); mpz_init(&a1); mpz_init(&a2); mpz_init(&gcd);
  mpz_init(&x1); mpz_init(&x2);
  mpz_init(&y1); mpz_init(&y2);
  mpz_init(&Ms); mpz_init(&Mi);
  for(i=0;i<n;i++)
   for(j=0;j<n; j++)
      mpz_set_si(&Trf[i][j], 0);
  for(i=0;i<n;i++)
      mpz_set_si(&Trf[i][i], 1);
  for(step = 0; step < n && spos != cols; step++)
  {
    tester = FALSE;
    while(tester == FALSE && spos != cols)
    {
       i = step;
       while(i<n && mpz_cmp_si(&M[i][spos], 0) == 0)
          i++;
       if(i<n)
       {
         tester = TRUE;
         if(i != step)
         {
           v = M[i];
           M[i] = M[step];
           M[step] = v;
           v = Trf[i];
           Trf[i] = Trf[step];
           Trf[step] = v;
         }
       }
       else
        spos++;
    }
    for(i=step+1;i<n && spos != cols;i++)
    {
      if(mpz_cmp_si(&M[i][spos], 0) != 0)
      {
        mpz_gcdext(&gcd, &a1, &a2, &M[step][spos], &M[i][spos]);
       if(mpz_cmp_si(&a1, 0) == 0)
       {
         v = M[step];
         M[step] = M[i];
         M[i] = v;
         v = Trf[step];
         Trf[step] = Trf[i];
         Trf[i] = v;
         mpz_set(&f, &a1); mpz_set(&a1, &a2); mpz_set(&a2, &f);
       }
       mpz_div(&Ms, &M[step][spos], &gcd);
       mpz_div(&Mi, &M[i][spos], &gcd);
       mpz_set(&M[step][spos], &Ms);
       mpz_set(&M[i][spos], &Mi);
       for(j=0;j<n;j++)
       {
         mpz_mul(&x1, &Trf[step][j], &a1);
         mpz_mul(&x2, &Trf[i][j], &a2);
         mpz_mul(&y1, &Trf[i][j], &Ms);
         mpz_mul(&y2, &Trf[step][j], &Mi);
         mpz_add(&Trf[step][j], &x1, &x2);
         mpz_sub(&Trf[i][j], &y1, &y2);
/**************************
         x = Trf->array.SZ[step][j];
         y = Trf->array.SZ[i][j];
         Trf->array.SZ[step][j] = a1 * x + a2 * y;
         Trf->array.SZ[i][j] = Ms * y - Mi * x;
**************************/
       }
       for(j=spos + 1; j<cols; j++)
       {
         mpz_mul(&x1, &M[step][j], &a1);
         mpz_mul(&x2, &M[i][j], &a2);
         mpz_mul(&y1, &M[i][j], &Ms);
         mpz_mul(&y2, &M[step][j], &Mi);
         mpz_add(&M[step][j], &x1, &x2);
         mpz_sub(&M[i][j], &y1, &y2);
/****************************
         x = M->array.SZ[step][j];
         y = M->array.SZ[i][j];
         M->array.SZ[step][j] = a1 * x + a2 * y;
         M->array.SZ[i][j] = Ms * y - Mi * x;
**************************/
       }
       mpz_set(&M[step][spos], &gcd);
       mpz_set_si(&M[i][spos], 0);
      }
    }
    if(spos != cols)
      rang++;
  }
  mpz_set_si(&y2, 2);
  mpz_set_si(&y1, 1);
  for(step = 0; step < rang; step++)
  {
    for(i=0;i<cols && mpz_cmp_si(&M[step][i], 0) == 0;i++);
    if(i != cols && mpz_cmp_si(&M[step][i], 0) < 0)
    {
       for(j=i;j<cols;j++)
         mpz_neg(&M[step][j], &M[step][j]);
       for(j=0;j<rows;j++)
         mpz_neg(&Trf[step][j], &Trf[step][j]);
    }
  }
  for(step = 0; step <rang-1;step++)
  {
     for(i=step+1;i<rang;i++)
     {
       for(spos = 0; spos<cols && mpz_cmp_si(&M[i][spos], 0) == 0;spos++);
       if(spos != cols)
       {
          mpz_set_si(&f, 0);
          mpz_div(&a1, &M[i][spos], &y2);
          if(mpz_cmp(&a1, &M[step][spos]) < 0)
          {
             mpz_div(&f, &M[step][spos], &M[i][spos]);
             mpz_mul(&x1, &f, &M[i][spos]);
             mpz_sub(&x2, &M[step][spos], &x1);
             if(mpz_cmp(&a1, &x2) < 0)
               mpz_add(&f, &f, &y1);
          }
          mpz_neg(&a1, &a1);
          if(mpz_cmp(&a1, &M[step][spos]) > 0)
          {
             mpz_div(&f, &M[step][spos], &M[i][spos]);
             mpz_mul(&x1, &f, &M[i][spos]);
             mpz_sub(&x2, &M[step][spos], &x1);
             if(mpz_cmp(&a1, &x2) > 0)
               mpz_sub(&f, &f, &y1);
             mpz_mod(&x1, &M[i][spos], &y2);
             if(mpz_cmp_si(&y2, 0) == 0 && mpz_cmp(&a1, &x2) == 0)
               mpz_sub(&f, &f, &y1);
          }
          if(mpz_cmp_si(&f, 0) != 0)
          {
             for(j=spos;j<cols;j++)
             {
                mpz_mul(&x1, &f, &M[i][j]);
                mpz_sub(&M[step][j], &M[step][j], &x1);
             }
             for(j=0;j<rows;j++)
             {
                mpz_mul(&x1, &f, &Trf[i][j]);
                mpz_sub(&Trf[step][j], &Trf[step][j], &x1);
             }
          }
       }
     }
  }
  
  mpz_clear(&f); mpz_clear(&a1); mpz_clear(&a2); mpz_clear(&gcd);
  mpz_clear(&x1); mpz_clear(&x2);
  mpz_clear(&y1); mpz_clear(&y2);
  mpz_clear(&Ms); mpz_clear(&Mi);
  return(rang);
}





/**************************************************************************\
@---------------------------------------------------------------------------
@ int MP_hnf(M, rows, cols)
@ MP_INT **M;
@ int rows, cols;
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
MP_hnf (MP_INT **M, int rows, int cols)
{
  int i,j;
  int n;
  int step;
  MP_INT a1,a2,gcd, *v, f;
  int tester;
  int spos = 0;
  MP_INT x1,x2,y1,y2;
  MP_INT Mi, Ms;
  int rang = 0;

  n = rows;
  mpz_init(&f); mpz_init(&a1); mpz_init(&a2); mpz_init(&gcd);
  mpz_init(&x1); mpz_init(&x2);
  mpz_init(&y1); mpz_init(&y2);
  mpz_init(&Ms);
  mpz_init(&Mi);

  for(step = 0; step < n && spos != cols; step++)
  {
    tester = FALSE;
    while(tester == FALSE && spos != cols)
    {
       i = step;
       while(i<n && mpz_cmp_si(&M[i][spos], 0) == 0)
          i++;
       if(i<n)
       {
         tester = TRUE;
         if(i != step)
         {
           v = M[i];
           M[i] = M[step];
           M[step] = v;
         }
       }
       else
        spos++;
    }
    for(i=step+1;i<n && spos != cols;i++)
    {
      if(mpz_cmp_si(&M[i][spos], 0) != 0)
      {
        mpz_gcdext(&gcd, &a1, &a2, &M[step][spos], & M[i][spos]);
       if(mpz_cmp_si(&a1, 0) == 0)
       {
         v = M[step];
         M[step] = M[i];
         M[i] = v;
         mpz_set(&f, &a1);
         mpz_set(&a1, &a2);
         mpz_set(&a2, &f);
       }
       mpz_div(&Ms, &M[step][spos], &gcd);
       mpz_div(&Mi, &M[i][spos], &gcd);
       mpz_set(&M[step][spos], &Ms);
       mpz_set(&M[i][spos], &Mi);
       for(j=spos + 1; j<cols; j++)
       {
         mpz_mul(&x1, &M[step][j], &a1);
         mpz_mul(&x2, &M[i][j], &a2);
         mpz_mul(&y1, &M[i][j], &Ms);
         mpz_mul(&y2, &M[step][j], &Mi);
         mpz_add(&M[step][j], &x1, &x2);
         mpz_sub(&M[i][j], &y1, &y2);
/****************************
         x = M->array.SZ[step][j];
         y = M->array.SZ[i][j];
         M->array.SZ[step][j] = a1 * x + a2 * y;
         M->array.SZ[i][j] = Ms * y - Mi * x;
**************************/
       }
       mpz_set(&M[step][spos], &gcd);
       mpz_set_si(&M[i][spos], 0);
      }
    }
    if(spos != cols)
      rang++;
  }
  mpz_set_si(&y2, 2);
  mpz_set_si(&y1, 1);
  for(step = 0; step < rang; step++)
  {
    for(i=0;i<cols && mpz_cmp_si(&M[step][i], 0) == 0;i++);
    if(i != cols && mpz_cmp_si(&M[step][i], 0) < 0)
    {
       for(j=i;j<cols;j++)
         mpz_neg(&M[step][j], &M[step][j]);
    }
  }
  for(step = 0; step <rang-1;step++)
  {
     for(i=step+1;i<rang;i++)
     {
       for(spos = 0; spos<cols && mpz_cmp_si(&M[i][spos], 0) == 0;spos++);
       if(spos != cols)
       {

          mpz_div(&f,&M[step][spos],&M[i][spos]);

          if (mpz_cmp_si(&M[step][spos],0) < 0){
             mpz_mod(&x1,&M[step][spos],&M[i][spos]);
             if (mpz_cmp_si(&x1,0) != 0)
                mpz_sub_ui(&f,&f,1);
          }

          /*
          mpz_set_si(&f, 0);
          mpz_div(&a1, &M[i][spos], &y2);
          if(mpz_cmp(&a1, &M[step][spos]) < 0)
          {
             mpz_div(&f, &M[step][spos], &M[i][spos]);
             mpz_mul(&x1, &f, &M[i][spos]);
             mpz_sub(&x2, &M[step][spos], &x1);
             if(mpz_cmp(&a1, &x2) < 0)
               mpz_add(&f, &f, &y1);
          }
          mpz_neg(&a1, &a1);
          if(mpz_cmp(&a1, &M[step][spos]) >= 0)
          {
             mpz_div(&f, &M[step][spos], &M[i][spos]);
             mpz_mul(&x1, &f, &M[i][spos]);
             mpz_sub(&x2, &M[step][spos], &x1);
             if(mpz_cmp(&a1, &x2) > 0)
               mpz_sub(&f, &f, &y1);
             mpz_mod(&x1, &M[i][spos], &y2);
             if(mpz_cmp_si(&y2, 0) == 0 && mpz_cmp(&a1, &x2) == 0)
               mpz_sub(&f, &f, &y1);
          } */
          if(mpz_cmp_si(&f, 0) != 0)
          {
             for(j=spos;j<cols;j++)
             {
                mpz_mul(&x1, &f, &M[i][j]);
                mpz_sub(&M[step][j], &M[step][j], &x1);
             }
             i--;
          }
       }
     }
  }
  mpz_clear(&f); mpz_clear(&a1); mpz_clear(&a2); mpz_clear(&gcd);
  mpz_clear(&x1); mpz_clear(&x2);
  mpz_clear(&y1); mpz_clear(&y2);
  mpz_clear(&Ms);
  mpz_clear(&Mi);

  return(rang);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ int MP_hnf_simultaneous(M, rows, cols, B, Bcols)
@ MP_INT **M, **B;
@ int rows, cols, Bcols;
@
@---------------------------------------------------------------------------
\**************************************************************************/
int 
MP_hnf_simultaneous (MP_INT **M, int rows, int cols, MP_INT **B, int Bcols)
{
  int i,j;
  int n;
  int step;
  MP_INT a1,a2,gcd, *v, f;
  int tester;
  int spos = 0;
  MP_INT x1,x2,y1,y2;
  MP_INT Mi, Ms;
  int rang = 0;

  n = rows;
  mpz_init(&f); mpz_init(&a1); mpz_init(&a2); mpz_init(&gcd);
  mpz_init(&x1); mpz_init(&x2);
  mpz_init(&y1); mpz_init(&y2);
  mpz_init(&Ms); mpz_init(&Mi);
  for(step = 0; step < n && spos != cols; step++)
  {
    tester = FALSE;
    while(tester == FALSE && spos != cols)
    {
       i = step;
       while(i<n && mpz_cmp_si(&M[i][spos], 0) == 0)
          i++;
       if(i<n)
       {
         tester = TRUE;
         if(i != step)
         {
           v = M[i];
           M[i] = M[step];
           M[step] = v;
           v = B[i];
           B[i] = B[step];
           B[step] = v;
         }
       }
       else
        spos++;
    }
    for(i=step+1;i<n && spos != cols;i++)
    {
      if(mpz_cmp_si(&M[i][spos], 0) != 0)
      {
        mpz_gcdext(&gcd, &a1, &a2, &M[step][spos], & M[i][spos]);
       if(mpz_cmp_si(&a1, 0) == 0)
       {
         v = M[step];
         M[step] = M[i];
         M[i] = v;
         v = B[step];
         B[step] = B[i];
         B[i] = v;
         mpz_set(&f, &a1); mpz_set(&a1, &a2); mpz_set(&a2, &f);
       }
       mpz_div(&Ms, &M[step][spos], &gcd);
       mpz_div(&Mi, &M[i][spos], &gcd);
       mpz_set(&M[step][spos], &Ms);
       mpz_set(&M[i][spos], &Mi);
       for(j=0;j<Bcols;j++)
       {
         mpz_mul(&x1, &B[step][j], &a1);
         mpz_mul(&x2, &B[i][j], &a2);
         mpz_mul(&y1, &B[i][j], &Ms);
         mpz_mul(&y2, &B[step][j], &Mi);
         mpz_add(&B[step][j], &x1, &x2);
         mpz_sub(&B[i][j], &y1, &y2);
/**************************
         x = Trf->array.SZ[step][j];
         y = Trf->array.SZ[i][j];
         Trf->array.SZ[step][j] = a1 * x + a2 * y;
         Trf->array.SZ[i][j] = Ms * y - Mi * x;
**************************/
       }
       for(j=spos + 1; j<cols; j++)
       {
         mpz_mul(&x1, &M[step][j], &a1);
         mpz_mul(&x2, &M[i][j], &a2);
         mpz_mul(&y1, &M[i][j], &Ms);
         mpz_mul(&y2, &M[step][j], &Mi);
         mpz_add(&M[step][j], &x1, &x2);
         mpz_sub(&M[i][j], &y1, &y2);
/****************************
         x = M->array.SZ[step][j];
         y = M->array.SZ[i][j];
         M->array.SZ[step][j] = a1 * x + a2 * y;
         M->array.SZ[i][j] = Ms * y - Mi * x;
**************************/
       }
       mpz_set(&M[step][spos], &gcd);
       mpz_set_si(&M[i][spos], 0);
      }
    }
    if(spos != cols)
      rang++;
  }
  mpz_set_si(&y2, 2);
  mpz_set_si(&y1, 1);
  for(step = 0; step < rang; step++)
  {
    for(i=0;i<cols && mpz_cmp_si(&M[step][i], 0) == 0;i++);
    if(i != cols && mpz_cmp_si(&M[step][i], 0) < 0)
    {
       for(j=i;j<cols;j++)
         mpz_neg(&M[step][j], &M[step][j]);
       for(j=0;j<Bcols;j++)
         mpz_neg(&B[step][j], &B[step][j]);
    }
  }
  for(step = 0; step <rang-1;step++)
  {
     for(i=step+1;i<rang;i++)
     {
       for(spos = 0; spos<cols && mpz_cmp_si(&M[i][spos], 0) == 0;spos++);
       if(spos != cols)
       {
          mpz_set_si(&f, 0);
          mpz_div(&a1, &M[i][spos], &y2);
          if(mpz_cmp(&a1, &M[step][spos]) < 0)
          {
             mpz_div(&f, &M[step][spos], &M[i][spos]);
             mpz_mul(&x1, &f, &M[i][spos]);
             mpz_sub(&x2, &M[step][spos], &x1);
             if(mpz_cmp(&a1, &x2) < 0)
               mpz_add(&f, &f, &y1);
          }
          mpz_neg(&a1, &a1);
          if(mpz_cmp(&a1, &M[step][spos]) > 0)
          {
             mpz_div(&f, &M[step][spos], &M[i][spos]);
             mpz_mul(&x1, &f, &M[i][spos]);
             mpz_sub(&x2, &M[step][spos], &x1);
             if(mpz_cmp(&a1, &x2) > 0)
               mpz_sub(&f, &f, &y1);
             mpz_mod(&x1, &M[i][spos], &y2);
             if(mpz_cmp_si(&y2, 0) == 0 && mpz_cmp(&a1, &x2) == 0)
               mpz_sub(&f, &f, &y1);
          }
          if(mpz_cmp_si(&f, 0) != 0)
          {
             for(j=spos;j<cols;j++)
             {
                mpz_mul(&x1, &f, &M[i][j]);
                mpz_sub(&M[step][j], &M[step][j], &x1);
             }
             for(j=0;j<Bcols;j++)
             {
                mpz_mul(&x1, &f, &B[i][j]);
                mpz_sub(&B[step][j], &B[step][j], &x1);
             }
          }
       }
     }
  }
  
  mpz_clear(&f); mpz_clear(&a1); mpz_clear(&a2); mpz_clear(&gcd);
  mpz_clear(&x1); mpz_clear(&x2);
  mpz_clear(&y1); mpz_clear(&y2);
  mpz_clear(&Ms); mpz_clear(&Mi);
  return(rang);
}
