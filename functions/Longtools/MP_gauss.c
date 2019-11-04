#include "typedef.h"
#include "gmp.h"
/* #include "gmp-impl.h" */
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: MP_gauss.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/************************************************************************\
@-------------------------------------------------------------------------
@ int MP_trf_gauss(M, Trf, rows, cols)
@ MP_INT **M, **Trf;
@ int rows, cols;
@ The Matrix M is transformed to a matrix M such that
@ M_new = Trf * M_old is Gauss-reduced, with Trf integral with
@ determinant +-1.
@ CAUTION: The entries of the matrix M are changed.
@-------------------------------------------------------------------------
\************************************************************************/

int 
MP_trf_gauss (MP_INT **M, MP_INT **Trf, int rows, int cols)
{
  int i,j;
  int n;
  int step;
  MP_INT a1,a2,gcd, *v, f;
  int tester;
  int spos = 0;
  MP_INT x1,x2,y1,y2;
  MP_INT Mi, Ms;

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
  for(step = 0; step < n; step++)
  {
    tester = FALSE;
    while(tester == FALSE)
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
       if(spos == cols)
       {
         mpz_clear(&f); mpz_clear(&a1); mpz_clear(&a2); mpz_clear(&gcd);
         mpz_clear(&x1); mpz_clear(&x2);
         mpz_clear(&y1); mpz_clear(&y2);
         mpz_clear(&Ms); mpz_clear(&Mi);
         return(step);
       }
    }
    for(i=step+1;i<n;i++)
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
  }
  mpz_clear(&f); mpz_clear(&a1); mpz_clear(&a2); mpz_clear(&gcd);
  mpz_clear(&x1); mpz_clear(&x2);
  mpz_clear(&y1); mpz_clear(&y2);
  mpz_clear(&Ms); mpz_clear(&Mi);
  return(n);
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ int MP_row_gauss(M, rows, cols)
@ MP_INT **M;
@ int rows, cols;
@
@ The same as MP_trf_gauss without calulating the transformation
@ matrix.
@---------------------------------------------------------------------------
@
\**************************************************************************/

int 
MP_row_gauss (MP_INT **M, int rows, int cols)
{
  int i,j;
  int n;
  int step;
  MP_INT a1,a2,gcd, *v, f;
  int tester;
  int spos = 0;
  MP_INT x1,x2,y1,y2;
  MP_INT Mi, Ms;

  n = rows;
  mpz_init(&f); mpz_init(&a1); mpz_init(&a2); mpz_init(&gcd);
  mpz_init(&x1); mpz_init(&x2);
  mpz_init(&y1); mpz_init(&y2);
  mpz_init(&Ms); mpz_init(&Mi);
  for(step = 0; step < n; step++)
  {
    tester = FALSE;
    while(tester == FALSE)
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
       if(spos == cols)
       {
         mpz_clear(&f); mpz_clear(&a1); mpz_clear(&a2); mpz_clear(&gcd);
         mpz_clear(&x1); mpz_clear(&x2);
         mpz_clear(&y1); mpz_clear(&y2);
         mpz_clear(&Ms); mpz_clear(&Mi);
         return(step);
       }
    }
    for(i=step+1;i<n;i++)
    {
      if(mpz_cmp_si(&M[i][spos], 0) != 0)
      {
        mpz_gcdext(&gcd, &a1, &a2, &M[step][spos], & M[i][spos]);
       if(mpz_cmp_si(&a1, 0) == 0)
       {
         v = M[step];
         M[step] = M[i];
         M[i] = v;
         mpz_set(&f, &a1); mpz_set(&a1, &a2); mpz_set(&a2, &f);
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
  }
  mpz_clear(&f); mpz_clear(&a1); mpz_clear(&a2); mpz_clear(&gcd);
  mpz_clear(&x1); mpz_clear(&x2);
  mpz_clear(&y1); mpz_clear(&y2);
  mpz_clear(&Ms); mpz_clear(&Mi);
  return(n);
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ int MP_row_gauss_simultaneous(M, rows, cols, B, Bcols)
@ MP_INT **M, **B;
@ int rows, cols, Bcols;
@
@ applies an integral gaussian algorithm (on the rows) to to 2-dimensional
@ array 'M' and does the simultaneous row operations on B.
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
MP_row_gauss_simultaneous (MP_INT **M, int rows, int cols, MP_INT **B, int Bcols)
{
  int i,j;
  int n;
  int step;
  MP_INT a1,a2,gcd, *v, f;
  int tester;
  int spos = 0;
  MP_INT x1,x2,y1,y2;
  MP_INT Mi, Ms;

  n = rows;
  mpz_init(&f); mpz_init(&a1); mpz_init(&a2); mpz_init(&gcd);
  mpz_init(&x1); mpz_init(&x2);
  mpz_init(&y1); mpz_init(&y2);
  mpz_init(&Ms); mpz_init(&Mi);
  for(step = 0; step < n; step++)
  {
    tester = FALSE;
    while(tester == FALSE)
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
       if(spos == cols)
       {
         mpz_clear(&f); mpz_clear(&a1); mpz_clear(&a2); mpz_clear(&gcd);
         mpz_clear(&x1); mpz_clear(&x2);
         mpz_clear(&y1); mpz_clear(&y2);
         mpz_clear(&Ms); mpz_clear(&Mi);
         return(step);
       }
    }
    for(i=step+1;i<n;i++)
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
  }
  mpz_clear(&f); mpz_clear(&a1); mpz_clear(&a2); mpz_clear(&gcd);
  mpz_clear(&x1); mpz_clear(&x2);
  mpz_clear(&y1); mpz_clear(&y2);
  mpz_clear(&Ms); mpz_clear(&Mi);
  return(n);
}


/*************************************************************************
@
@------------------------------------------------------------------------
@
@ void MP_row_gauss_reverse(MP_INT **A,int rows,int cols,int option)
@
@ Performs a manhattan style gauss algorithm on the MP_INT matrix
@ M, which has to be gauss reduced by MP_row_gauss before.
@ The algorithm is not done complete for sake of speed!
@
@ MP_INT **A   : The matrix in question. It will be changed!
@ int rows     : the rows of A. It is suffisient to feed the rank in here.
@ int cols     : the number of columns of A.
@ int option   : this flag determines the behaviour of the function:
@                for option == 0 it will only perform Z elementary
@                tranformations, for option == 1 it will divide every
@                row by the gcd of it's entries (so do it Q style)!
@------------------------------------------------------------------------
@
*************************************************************************/
void MP_row_gauss_reverse(MP_INT **A,int rows,int cols,int option)
{
   int i,
       j,
       k,
       actcol;

   MP_INT x,
          y;

   mpz_init(&x);
   mpz_init(&y);

   for (k=rows-1;k>=0;k--){

      /* find the first non zero entry in the k-th row of A */
      for (actcol=0;actcol<cols && (mpz_cmp_si(A[k]+actcol,0) ==0);actcol++);

      /* make the first entry in this row positive for nice purposes */
      if (mpz_cmp_si(A[k]+actcol,0)<0)
        for (j=actcol;j<cols;j++) mpz_neg(A[k]+j,A[k]+j);

      /* bare in mind: option == 1 means we are doing gauss over Q, so
         divide the row by the gcd of it's entries */
      if (option == 1){
         mpz_set(&x,A[k]+actcol);
         for (i=actcol+1;i<cols;i++) mpz_gcd(&x,&x,A[k]+i);
         mpz_abs(&x,&x);
         if (mpz_cmp_si(&x,1) != 0)
            for (i=actcol;i<cols;i++) mpz_div(A[k]+i,A[k]+i,&x);
      }

      for (i=0;i<k;i++){
         mpz_div(&x,A[i]+actcol,A[k]+actcol);
         for (j=actcol;j<cols;j++){
            mpz_mul(&y,&x,A[k]+j);
            mpz_sub(A[i]+j,A[i]+j,&y);
         }
      }
   }

   mpz_clear(&x);
   mpz_clear(&y);

   return;
}

