#include"typedef.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  reduction_sort.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ void reduction_sort(G,T,n)
@ int **G, **T, n;
@
@ make simultaneous row and colum permutations on the twodimensional
@ array G of size n x n such that G[i][i] <= G[j][j] for 0<=i<=j<n.
@ Simultaneously the row operations are applied to T, that
@ has to be an array of the same size or NULL
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
reduction_sort (int **G, int **T, int n)
{
  int i,j;
  int minpos, min;
  int *tmp;
  int merk;
  for(i=0;i<n;i++)
  {
    min = G[i][i];
    minpos = i;
    for(j=i+1;j<n;j++)
    {
     if(G[j][j] < min)
     { min = G[j][j]; minpos = j;}
    }
    if(minpos != i)
    {
       tmp = G[i]; G[i] = G[minpos]; G[minpos] = tmp;
       if(T != NULL)
       { tmp = T[i]; T[i] = T[minpos]; T[minpos] = tmp;}
       for(j=0;j<n;j++)
       { merk = G[j][i]; G[j][i] = G[j][minpos]; G[j][minpos] = merk;}
    }
  }
}
