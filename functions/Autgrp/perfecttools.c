#include "typedef.h"
#include "utils.h"
#include "types.h"

/*******************************************************************\
| calculate the partial image of the matrices in 'perpbase'
| under the partial automorphism and check, whether it is
| equal to the IxI-submatrix (of the first I columns and rows)
| of one of the matrices in 'perp'
\*******************************************************************/


int normal_aut_test (int *x, int I, veclist V)
{
  int i,j,k,l;
  int o,u,m;
  int found, tester;
  int xl, xj;

  for(i=0;i<perpdim;i++)
  {
    for(l=0;l<=I;l++)
    {
      if(x[l] < 0) xl = -x[l]; else xl = x[l];
      /***************************************************************\
      | calculate the l-th row of perpprod[i]
      \***************************************************************/
      for(j=0;j<V.dim;j++)
      {
         perpvec[j] = 0;
         for(k=0;k<V.dim;k++)
           perpvec[j] += V.v[xl][k] * perpbase[i][k][j];
      }
      if(x[l] < 0)
      {
          for(j=0;j<V.dim;j++)
            perpvec[j] = -perpvec[j];
      }
      for(j=0;j<=l;j++)
      {
        perpprod[i][l][j] = 0;
        if(x[j] < 0) xj = -x[j]; else xj = x[j];
        for(k=0;k<V.dim;k++)
          perpprod[i][l][j] += perpvec[k] * V.v[xj][k];
        if(x[j] < 0)
          perpprod[i][l][j] = -perpprod[i][l][j];
      }
    }
    /******************************************************************\
    | Check, whether perprod[i] is a submatrix of a matrix in perp
    | The matrices in perp are assumed to be sorted.
    \******************************************************************/
    found = FALSE;
    u = 0; o = perp_no-1;
    while(found == FALSE && o >=u )
    {
      m = (o+u)/2;
      tester = 0;
      for(j=0;j<=I && tester == 0;j++)
      {
          for(k=0;k<=j && tester == 0;k++)
          {
            if(perp[m][j][k] != perpprod[i][j][k])
            {
              if(perp[m][j][k] > perpprod[i][j][k])
                 tester = 1;
              else tester = -1;
            }
          }
      }
      if(tester == 0)
         found = TRUE;
      else
      {
         if(tester == 1) o = m-1;
         else u = m+1;
      }
    }
    if(found == FALSE)
    {
      return(FALSE);
    }
  }
  return(TRUE);
}


void mach_perp_matrices (fpstruct fp, matrix_TYP **P, matrix_TYP **Pbase, int n)
{
        int i,j,k;

        perp = (int ***)xmalloc(perp_no *sizeof(int **));
        for(i=0;i<perp_no;i++)
        {
         perp[i] = (int **)xmalloc(n *sizeof(int *));
         for(j=0;j<n;j++)
         {
          perp[i][j] = (int *)xmalloc(n *sizeof(int));
         }
         for(j=0;j<n;j++)
           for(k=0;k<n;k++)
            perp[i][j][k] = P[i]->array.SZ[fp.per[j]][fp.per[k]];
        }
        perpbase = (int ***)xmalloc(perpdim *sizeof(int **));
        for(i=0;i<perpdim;i++)
        {
         perpbase[i] = (int **)xmalloc(n *sizeof(int *));
         for(j=0;j<n;j++)
         {
          perpbase[i][j] = (int *)xmalloc(n *sizeof(int));
         }
         for(j=0;j<n;j++)
           for(k=0;k<n;k++)
            perpbase[i][j][k] = Pbase[i]->array.SZ[j][k];
        }
        perpprod = (int ***)xmalloc(perpdim *sizeof(int **));
        for(i=0;i<perpdim;i++)
        {
          perpprod[i] = (int **)xmalloc(n *sizeof(int *));
          for(j=0;j<n;j++)
          {
            perpprod[i][j] = (int *)xmalloc((j+1) *sizeof(int));
          }
        }
        perpvec = (int *)xmalloc(n *sizeof(int));
}

void free_perp_matrices (int n)
{
   int i,j;
   for(i=0;i<perp_no;i++)
   {
      for(j=0;j<n;j++)
        free(perp[i][j]);
      free(perp[i]);
   }
   free(perp);
   for(i=0;i<perpdim;i++)
   {
      for(j=0;j<n;j++)
        free(perpbase[i][j]);
      free(perpbase[i]);
   }
   free(perpbase);
   for(i=0;i<perpdim;i++)
   {
      for(j=0;j<n;j++)
        free(perpprod[i][j]);
      free(perpprod[i]);
   }
   free(perpprod);
   free(perpvec);
}
