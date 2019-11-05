#include "typedef.h"
#include "matrix.h"
#include "tools.h"

static int n;

/************************************************************************\
| The Matrix M is transformed to a matrix M' such that
| M' = Trf * M is Gauss-reduced, with Trf integral with determinante +-1.
\************************************************************************/

int 
Trf_gauss (matrix_TYP *M, matrix_TYP *Trf)
{
  int i,j;
  int step;
  int a1,a2,gcd;
  int tester, *v;
  int spos = 0;
  int x,y;
  int Mi, Ms;

  if(Trf->cols != M->rows)
  {
       printf("Error in Trf_gauss: matrices M and Trf are not compatible\n");
       exit(3);
  }
  if(Trf->cols != Trf->rows)
  {
       printf("Error in Trf_gauss: Trf must be a square matrix\n");
       exit(3);
  }
  n = Trf->rows;
  for(i=0;i<n;i++)
   for(j=0;j<n; j++)
      Trf->array.SZ[i][j] = 0;
  for(i=0;i<n;i++)
      Trf->array.SZ[i][i] = 1;
  for(step = 0; step < n; step++)
  {
    tester = FALSE;
    while(tester == FALSE)
    {
       i = step;
       while(i<n && M->array.SZ[i][spos] == 0)
          i++;
       if(i<n)
       {
         tester = TRUE;
         if(i != step)
         {
           v = M->array.SZ[i];
           M->array.SZ[i] = M->array.SZ[step];
           M->array.SZ[step] = v;
           v = Trf->array.SZ[i];
           Trf->array.SZ[i] = Trf->array.SZ[step];
           Trf->array.SZ[step] = v;
         }
       }
       else
        spos++;
       if(spos == M->cols)
         return(step);
    }
    for(i=step+1;i<n;i++)
    {
      if(M->array.SZ[i][spos] != 0)
      {
        gcd_darstell(M->array.SZ[step][spos], M->array.SZ[i][spos], &a1, &a2, &gcd);
       if(a1 == 0)
       {
         v = M->array.SZ[step];
         M->array.SZ[step] = M->array.SZ[i];
         M->array.SZ[i] = v;
         v = Trf->array.SZ[step];
         Trf->array.SZ[step] = Trf->array.SZ[i];
         Trf->array.SZ[i] = v;
         j = a1; a1 = a2; a2 = j;
       }
       Ms = M->array.SZ[step][spos] / gcd;
       Mi = M->array.SZ[i][spos] / gcd;
       M->array.SZ[step][spos] /= gcd;
       M->array.SZ[i][spos] /= gcd;
       for(j=0;j<n;j++)
       {
         x = Trf->array.SZ[step][j];
         y = Trf->array.SZ[i][j];
         Trf->array.SZ[step][j] = a1 * x + a2 * y;
         Trf->array.SZ[i][j] = Ms * y - Mi * x;
       }
       for(j=spos + 1; j<M->cols; j++)
       {
         x = M->array.SZ[step][j];
         y = M->array.SZ[i][j];
         M->array.SZ[step][j] = a1 * x + a2 * y;
         M->array.SZ[i][j] = Ms * y - Mi * x;
       }
       M->array.SZ[step][spos] = gcd;
       M->array.SZ[i][spos] = 0;
      }
    }
  }
  return(n);
}

/************************************************************************\
| solve_mat(M) calculates an Matrix X with MX^{tr} = 0, such that
| the rows of X are a Z-basis of the solution space.
\************************************************************************/
matrix_TYP *
solve_mat (matrix_TYP *M)
{
   matrix_TYP *M1, *M1t, *Trf, *X, *erg;
   int i,rang, n;

   M1 = copy_mat(M);
   row_gauss(M1);
   M1t = tr_pose(M1);
   n = M1t->rows;
   Trf = init_mat(n, n, "");
   rang = Trf_gauss(M1t, Trf);
   X = init_mat(n-rang, n, "");
   for(i=0;i<n-rang;i++)
     X->array.SZ[i][rang+i] = 1;
   erg = mat_mul(X, Trf);
   free_mat(X); free_mat(Trf);
   free_mat(M1); free_mat(M1t);
   return(erg);
}
