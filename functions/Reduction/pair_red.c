#include"typedef.h"
#include"tools.h"
#include"reduction.h"

/************************************************************************** \
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: pair_red.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


static void two_reduce(G, T, i, j, n)
int **G, **T, i, j, n;
{
   int k, f;
  
   if(G[i][j] > 0)
     f = (G[i][j] + G[i][i]/2)/G[i][i];
   else
     f = (G[i][j] - G[i][i]/2)/G[i][i];
   
   for(k=0;k<n;k++)
       G[j][k] -=  f * G[i][k];
   G[j][j] -= f * G[j][i];
   for(k=0;k<n;k++)
   {
     if(k != j)
     {
       G[k][j] = G[j][k];
     }
   }
   for(k=0;k<n && T != NULL;k++)
       T[j][k] -=  f * T[i][k];
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ void pr_red(G, T, n)
@ int **G, **T, n;
@
@  applies a pair_reduction to the 2-dimensional array G of size n x n
@  and makes the simutaneous row operations on the array T of the same size
@ The function can be called with T = NULL
@---------------------------------------------------------------------------
@
\**************************************************************************/
void pr_red(G, T, n)
int **G, **T, n;
{
   int i,j,k, a;
   int reduced, red2;
   int merk, s, p1, p2;

   reduced = FALSE;
   reduction_sort(G,T,n);
   while(reduced == FALSE)
   {
     reduced = TRUE;
     for(i=0;i<n;i++)
      for(j=0;j<i;j++)
      {
         red2 = FALSE;
         while(red2 == FALSE && (G[i][i] != 0 || G[j][j] != 0))
         {
           if(G[i][i] > 0 && G[i][i] <=  G[j][j])
           {
             a = 2*G[i][j];
             if(abs(a) > G[i][i])
             {
               two_reduce(G, T, i, j, n);
               reduced = FALSE;
             }
             else
              red2 = TRUE;
           }
           if(G[j][j] > 0 && G[j][j] <= G[i][i])
           {
             a = 2*G[i][j];
             if(a > G[j][j] || -a > G[j][j])
             {
               two_reduce(G, T, j, i, n);
               reduced = FALSE;
             }
             else
              red2 = TRUE;
           }
           if(G[i][i] < 0 || G[j][j] < 0)
             return;
         }
      }
      reduction_sort(G, T, n);
      for(i=0;i<n;i++)
        for(j=i+1;j<n;j++)
        {
          if(G[i][i] <= G[j][j])
            { p1 = i; p2 = j;}
          else
            { p1 = j; p2 = i;}
          if(2*abs(G[p1][p2]) == G[p1][p1])
          {
            s = signum(G[p1][p2]);
            for(k=0; k<n;k++)
            {
              if(k != p2)
                G[p2][k] -= s * G[p1][k];
              if(T != NULL)
              T[p2][k] -= s * T[p1][k];
            }
            for(k=0;k<n;k++)
            {
              G[k][p2] = G[p2][k];
            }
            for(k=0;k<n;k++)
            {
              if(k != p1 && k != p2)
              {
                red2 = FALSE;
                while(red2 == FALSE && G[k][k] > 0 && G[p2][p2] > 0)
                {
                  a = 2*abs(G[k][p2]);
                  if(a > G[k][k] || a > G[p2][p2])
                  {
                     if(a > G[k][k] && G[k][k] < G[p2][p2])
                      two_reduce(G, T, k, p2, n);
                     else
                      two_reduce(G, T, p2, k, n);
                     reduced = FALSE;
                  }
                  else
                   red2 = TRUE;
                }
                if(G[k][k] < 0 || G[p2][p2] < 0)
                  return;
              }
            }
          }
          if(G[i][i] < 0 || G[j][j] < 0)
             return;
        }
   }
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *pair_red(Gram, Tr)
@ matrix_TYP *Gram, *Tr;
@
@ calculates matrices A and Tr such that A = Tr * Gram * Tr^{tr}
@ is pair_reduced.
@ The function can be called with Tr = NULL
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *pair_red(Gram, Tr)
matrix_TYP *Gram, *Tr;
{
  int i,j;
  int **G, **T;
  matrix_TYP *Gerg;

  extern matrix_TYP *copy_mat();

  Gerg = copy_mat(Gram);
  G = Gerg->array.SZ;
  if(Tr != NULL)
  {
    T = Tr->array.SZ;
    for(i=0;i<Gram->cols;i++)
    {
     for(j=0;j<Gram->cols;j++)
       T[i][j] = 0;
     T[i][i] = 1;
    }
  }
  else
    T = NULL;
  pr_red(G, T, Gram->cols);
  return(Gerg);
}
