#include "typedef.h"
#include "utils.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: dsylv.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/* changed tilman 2/09/97 because this is far higher than the precision
   of modern machines (and causes trouble)
static  double eps = 0.000001; */
static  double eps = 0.00000000001;

static double **
mat_to_double (matrix_TYP *M)
{
  int i,j;
  double **D;
  D = (double **)xmalloc(M->rows *sizeof(double *));
  for(i=0;i<M->rows;i++)
  {
    D[i] = (double *)xmalloc((i+1) *sizeof(double));
  }
  for(i=0;i<M->rows;i++)
    for(j=0;j<=i;j++)
        D[i][j] = ((double) M->array.SZ[i][j]);
  return(D);
}


static int 
pivot (double **A, int step, int n)
{
  int i,j, found;
  double max;
  int ipos, jpos;
  double tmp;
  
/*****************************************************************************\
|  search for a nonzero diagonal entry
\*****************************************************************************/
  max = fabs(A[step][step]);
  if (max < eps){
     max = A[step][step] = 0.0;
  }
  ipos = step;
  for(i=step+1;i<n;i++)
  {

    /* inserted tilman 02/09/97, accuracy errors may occur intermediate */
    if((fabs(A[i][i]) < eps))
       A[i][i] = 0.0;
    if((fabs(A[i][i]) > max))
    {
       max = fabs(A[i][i]);
       ipos = i;
    }
  }
/*****************************************************************************\
|  if there is a nonzero diagonal entry, swap colum 'step' und colmn 'ipos'
|  and swap row 'step' and row 'ipos', where A[ipos][ipos] has the maximal
|  absulute value.
\*****************************************************************************/
  if(max != 0.0)
  {
     if (ipos == step)
       return(1);
     tmp = A[step][step];
     A[step][step] = A[ipos][ipos];
     A[ipos][ipos] = tmp;
     for(i=step+1;i<ipos;i++)
         { tmp = A[i][step]; A[i][step] = A[ipos][i]; A[ipos][i] = tmp;}
     for(i= ipos+1;i<n;i++)
         { tmp = A[i][step]; A[i][step] = A[i][ipos]; A[i][ipos] = tmp;}
     return(1);
  }

/*****************************************************************************\
| if all diagonal entries are zero, find an nonzero entry A[ipos][jpos]
| if all entries are zero return(0)
\*****************************************************************************/
  ipos = step+1;
  jpos = step;
  found = FALSE;
  for(i=step+1;i<n && found == FALSE;i++)
    for(j=step;j<i && found == FALSE;j++)
    {
       if(A[i][j] != 0.0)
       {ipos = i; jpos = j; found = TRUE;}
    }
  if(found == FALSE)
    return(0);

/*****************************************************************************\
| swap colum 'step' and column 'jpos' and swap row 'step' and row 'jpos'
| Then A[ipos][step] is nonzero.
\*****************************************************************************/
     for(i=step+1;i<jpos;i++)
         { tmp = A[i][step]; A[i][step] = A[jpos][i]; A[jpos][i] = tmp;}
     for(i= jpos+1;i<n;i++)
         { tmp = A[i][step]; A[i][step] = A[i][jpos]; A[i][jpos] = tmp;}
/*****************************************************************************\
| add column 'ipos' to column 'step' and add row 'ipos' to row 'step'
| Then A[step][step] is nonzero
\*****************************************************************************/
   A[step][step] = A[ipos][step] * 2.0;
   for(i=step+1;i<ipos;i++)
     A[i][step] += A[ipos][i];
   for(i=ipos+1;i<n;i++)
     A[i][step] += A[i][ipos];
  return(1);
}


static void 
col_clear (double **A, int step, int n)
{
   int i,j;
   for(i=step+1; i<n; i++)
   {
     /* changed this to a numerically more stable version:
     if(fabs(A[i][step]) > eps)
     {
        double f = A[i][step]/A[step][step];
        for(j=step+1;j<=i;j++)
           A[i][j] -= A[j][step] * f;
        for(j=i+1;j<n;j++)
           A[j][i] -= A[j][step] * f;
        A[i][step] = 0.0;
     }
     else
     {
        A[i][step] = 0.0;
     } */

     for(j=step+1;j<=i;j++)
        A[i][j] -= (A[j][step] * A[i][step])/A[step][step];
     for(j=i+1;j<n;j++)
        A[j][i] -= (A[j][step] * A[i][step])/A[step][step];
     A[i][step] = 0.0;

   }
}


static int 
diag_definite_test (double **D, int n)
{
  int i;
  int o,u,z;
  
  o = u = z = 0;
  for(i=0;i<n;i++)
  {
     if(fabs(D[i][i]) > eps)
     {
       if(D[i][i] > eps)
         o = 1;
       else
         u = 1;
     }
     else 
       z = 1;
  }
  if(u == 0)
  {
    if(o == 0)
       return(0);
    if(z == 0)
       return(2);
    return(1);
  }
  if(o == 1)
    return(-3);
  if(z == 0)
    return(-2);
  return(-1);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *dsylv(M)
@ matrix_TYP *M;
@ 
@ The return of 'dsylv' is a diagonal matrix D with diaogonal elements
@ D[i][i] in {-1,1,0}.
@ The function diagonalizes the symmetric matrix 'M' by simultanous row
@ and col-operations.
@ If the entry in the diagonalised matrix is positiv it is replaced
@ by 1, if negativ by -1, if zero by 0.
@ This is done with floating-point arithmetic.
@ So rounding errows may occur.
@ Therefore a diaogonal entry it is assumed to be zero
@ if its absolute value is less then eps.
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
dsylv (matrix_TYP *M)
{
  int i,n,step, nonzero;
  double **D;
  matrix_TYP *A;
  int pos, neg; 

  D = mat_to_double(M);
  n = M->cols;
  nonzero = pivot(D, 0, n);
  for(step = 0;step<n-1 && nonzero == TRUE ;step++)
  {
       col_clear(D, step, n);

       if(step+1 != n)
         nonzero = pivot(D, step+1, n);
  }
  pos = neg = 0;
  for(i=0;i<n;i++)
  {
     if(fabs(D[i][i]) > eps)
     {
         if(D[i][i] > eps)
           pos++;
         else
           neg++;
     }
  }
  A = init_mat(n,n,"");
  A->flags.Symmetric = A->flags.Diagonal = TRUE;
  for(i=0;i<pos;i++)
     A->array.SZ[i][i] = 1;
  for(i=pos; i<pos+neg;i++)
     A->array.SZ[i][i] = -1;
  for(i=0;i<n;i++)
    free(D[i]);
  free(D);
  return(A);
}




/**************************************************************************\
@---------------------------------------------------------------------------
@ int definite_test(M)
@ matrix_TYP *M;
@ 
@ Does the same as 'dsolve' but testes during the Diagonalisation
@ if the matrix is indefinite.
@ The return is
@ 
@      0:  if M is zero
@      1:  if M is positiv semidefinite, but not positiv definite.
@      2:  if M is positiv definite.
@     -1:  if M is negativ semidefinite, but not negativ definite.
@     -2:  if M is negativ definite.
@     -3:  if M is indefinite.
@
@---------------------------------------------------------------------------
@
\**************************************************************************/

int 
definite_test (matrix_TYP *M)
{
  int i,n,step, nonzero;
  double **D;
  int test;

  
  D = mat_to_double(M);
  n = M->cols;
  test = diag_definite_test(D, n);  
  if(test == -3)
  {
    for(i=0;i<n;i++)
       free(D[i]);

    /* inserted the next line 15/1/97 tilman */
    free(D);

    return(-3);
  }
  nonzero = pivot(D, 0, n);
  for(step = 0;step<n-1 && nonzero == TRUE ;step++)
  {
       col_clear(D, step, n);
       test = diag_definite_test(D, n);  
       if(test == -3)
       {
         for(i=0;i<n;i++)
            free(D[i]);

         /* inserted the next line 15/1/97 tilman */
         free(D);

         return(-3);
       }
       if(step+1 != n)
         nonzero = pivot(D, step+1, n);
  }
  for(i=0;i<n;i++)
    free(D[i]);

  /* inserted the next line 15/1/97 tilman */
  free(D);

  return(test);
}
