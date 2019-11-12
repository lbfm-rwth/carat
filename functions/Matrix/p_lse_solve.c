#include "typedef.h"
#include "utils.h"
#include "matrix.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: p_lse_solve.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
/********************************************************************\
@  The function 'matrix_TYP **p_lse_solve(A, B, anz, p)'
@  calculates for given matrices A and B the solutions of 
@              A * X = B  (modulo p)
@  where A is a (n x m)- matrix and B a (n x r) matrix with integral
@  entries (treated as elements of the field with p elements).
@  The number of the returned matrices is given back by the pointer 'anz'.
@
@  If no solution exists a NULL-pointer is returned, and 'anz' is
@  set to 0.
@
@  If a solution of the inhomogenous equation exists, it is
@  returned via the (m x r)-matrix X[0] ( X is the matrix_TYP ** returned by
@  the functions.
@  The (1 x m)-matrices X[1],...,X[anz-1] form a basis of the solutions of
@  the transposed homogenous equation X * A^{tr} = 0.
@  
\********************************************************************/
matrix_TYP **
p_lse_solve (matrix_TYP *A, matrix_TYP *B, int *anz, int p)
{
  matrix_TYP **X;
  int **C, **D;
  int i,j,k,step, cpos, f;
  int m,n,r;

  int inhomo,
      found,
      rang = 0,   /* inserted tilman 21/07/97: if the matrix is zero mod p,
                     a different setting causes trouble */
      Xsize;

  int *tmp, *pos;
  int phalbe, mphalbe;
  
  phalbe = p/2;
  mphalbe = -phalbe;
  if(p == 2)
    mphalbe = 0;
  n = A->rows, m = A->cols;
  C = (int **)xmalloc(n *sizeof(int *));
  for(i=0;i<n;i++)
  {
    C[i] = (int *)xmalloc(m *sizeof(int));
    for(j=0;j<m;j++)
      C[i][j] = A->array.SZ[i][j];
  }
  if(B == NULL)
  {
    inhomo = FALSE;
    r = 0;
    D = 0;
  }
  else
  {
     inhomo = TRUE;
     r = B->cols;
     D = (int **)xmalloc(n *sizeof(int *));
     for(i=0;i<n;i++)
     {
       D[i] = (int *)xmalloc(r *sizeof(int));
       for(j=0;j<r;j++)
         D[i][j] = B->array.SZ[i][j];
     }
  }
  if(inhomo && n != B->rows)
  {
     printf("You tried to solve A * X = B, where A has %d rows and B has %d rows\n", A->rows, B->rows);
     printf("exit in 'p_lse_solve':\n");
     exit(3);
  }

  /*************************************************************************\
  |  reduce the entries in C and D modulo p
  \*************************************************************************/
   for(i=0;i<n;i++)
    for(j=0;j<m;j++)
    {
      C[i][j] %= p;
      if(C[i][j] > phalbe)
        C[i][j] -= p;
      if(C[i][j] < mphalbe)
        C[i][j] += p;
    }
   if(inhomo)
   {
     for(i=0;i<n;i++)
      for(j=0;j<r;j++)
      {
        D[i][j] %= p;
        if(D[i][j] > phalbe)
          D[i][j] -= p;
        if(D[i][j] < mphalbe)
          D[i][j] += p;
      }
   }

  /*************************************************************************\
  | Make a row-gauss-algorithm on C and do simultanous row operations on D.
  \*************************************************************************/
  cpos = 0;
  for(step = 0; step < n; step++)
  {
     /**************************************************************\
     | search first column 'cpos' in C with nonzero entry C[i][cpos]
     | and permute the rows 'step' and 'i'
     \**************************************************************/
     found = FALSE;
     while(cpos < m && found == FALSE)
     {
       for(i=step; i<n && C[i][cpos] == 0;i++);
       if(i == n)
        cpos++;
       else
         found = TRUE;
     }
     if(found == TRUE)
     {
       rang = step+1;
       if(i!= step)
       { 
         tmp = C[i]; C[i] = C[step]; C[step] = tmp; tmp = NULL;
         if(inhomo)
         {  tmp = D[i]; D[i] = D[step]; D[step] = tmp; tmp = NULL;}
       }
       /******************************************************************\
       | Clear the entries C[step+1][cpos],....,C[n][cpos]
       \******************************************************************/
       f = p_inv(C[step][cpos], p);
       if(f > phalbe)  f -= p;
       if(f < mphalbe)  f += p;
       if(f != 1)
       {
         C[step][cpos] = 1;
         for(i=cpos+1;i<m; i++)
         {
           C[step][i] = (C[step][i] * f)%p;
           if(C[step][i] > phalbe) C[step][i] -=p;
           if(C[step][i] < mphalbe) C[step][i] +=p;
         }
         C[step][cpos] = 1;
         if(inhomo)
         {
            for(i=0;i<r;i++)
            {
              D[step][i] = (D[step][i] * f)%p;
              if(D[step][i] > phalbe) D[step][i] -=p;
              if(D[step][i] < mphalbe) D[step][i] +=p;
            }
         }
       }
       for(i=step+1;i<n;i++)
       {
         if(C[i][cpos] != 0)
         {
            f = C[i][cpos];
            C[i][cpos] = 0;
            for(j=cpos+1; j<m;j++)
            {
             C[i][j] = ( C[i][j] - f * C[step][j])%p;
              if(C[i][j] > phalbe) C[i][j] -=p;
              if(C[i][j] < mphalbe) C[i][j] +=p;
            }
            if(inhomo)
            {
              for(j=0;j<r;j++)
              {
                D[i][j] = ( D[i][j] - f * D[step][j])%p;
                if(D[i][j] > phalbe) D[i][j] -=p;
                if(D[i][j] < mphalbe) D[i][j] +=p;
              }
            }
         }
       }
     }
  }

  /********************************************************************\
  | Now C is an upper triangular matrix and the number of nonzero rows
  | is given by the variable 'rang'.
  | If A[i][j] is the first non-zero entry in the i-th row (which has been
  | made equal to 1),
  | all the entries A[k][j] with k<i are made zero.
  \********************************************************************/
  for(i = rang-1; i >= 0; i--)
  {
    for(j=0; j<m && C[i][j] == 0; j++);
    cpos = j;
    for(j=0; j<i;j++)
    {
      f = C[j][cpos];
      for(k=cpos+1;k<m;k++)
      {
       C[j][k] = (C[j][k] - f*C[i][k])%p;
       if(C[j][k] > phalbe)  C[j][k] -= p;
       if(C[j][k] < mphalbe)  C[j][k] += p;
      }
      C[j][cpos] = 0;
      if(inhomo)
      {
        for(k=0;k<r;k++)
        {
         D[j][k] = (D[j][k] - f*D[i][k])%p;
         if(D[j][k] > phalbe)  D[j][k] -= p;
         if(D[j][k] < mphalbe)  D[j][k] += p;
        }
      }
    }
  }
  /********************************************************************\
  | Calculate a solution of the inhomogenous equation
  \********************************************************************/
  if(inhomo)
  {
     /********************************************************************\
     | Check, whether a solution of the inhomogenous equation exists.
     | If not return NULL-pointer
     \********************************************************************/
     found = FALSE;
     for(i=rang;i<n && found == FALSE;i++)
      for(j=0;j<r && found == FALSE; j++)
      {
        if(D[i][j] != 0) found = TRUE;
      }
     if(found == TRUE)
     {
       for(i=0;i<n;i++)
       {  free(D[i]); free(C[i]);}
       free(C); free(D);
       *anz = 0;
       return(NULL);
     }
  }
  Xsize = m-rang+1;
  *anz = Xsize;
  X = (matrix_TYP **)xmalloc(Xsize *sizeof(matrix_TYP *));

  /*******************************************************************\
  | pos[i] is the index j such that C[i][j] is the first nonzero
  | entry in the i-th row of C
  \*******************************************************************/
  pos = NULL;
  pos = (int *)xmalloc((rang+2) *sizeof(int));

  /* inserted tilman 11/4/97 */
  pos++; pos[-1] = -1;

  for(i=0;i<rang;i++)
  {
    for(j=0;j<m && C[i][j] == 0; j++);
    pos[i] = j;
  }
  pos[rang] = m;
  if(inhomo)
  {
    /****************************************************************\
    | Find a solution of the inhomogenous equation
    \****************************************************************/
    X[0] = init_mat(m, r, "");
    for(i=0;i<r;i++)
      for(j=0;j<rang;j++)
        X[0]->array.SZ[pos[j]][i] = D[j][i];
    for(i=0;i<n;i++)
     free(D[i]);
    free(D);
  }
  else
   X[0] = NULL;
  /****************************************************************\
  | Find the solutions of the homogenous equation
  \****************************************************************/
  step=1;

  /* changed 11/4/97 tilman from
  for(i=rang-1; i>= 0; i--) to */
  for(i=rang-1; i>= (-1); i--)
  {
    for(j=pos[i]+1; j< pos[i+1]; j++)
    {
       X[step] = init_mat(1, m, "");
       tmp = X[step]->array.SZ[0];
       step++;
       tmp[j] = -1;
       for(k=0;k<=i;k++)
          tmp[pos[k]] = C[k][j];
    }
  }
  for(i=0;i<n;i++)
   free(C[i]);
  free(C);

  /* inserted tilman 11/4/97 */
  free(--pos);
  return(X);
}
