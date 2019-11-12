#include "typedef.h"
#include "utils.h"
#include"matrix.h"
#include"longtools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: formspace.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP **formspace(B, Banz, sym_opt, fdim)
@ matrix_TYP **B;
@ int Banz, sym_opt, *fdim;
@
@  formspace calculates a basis of the lattice of integral matrices X with
@      B[i]^{tr} * X * B[i] = X    for 0 <= i < Banz
@  The number of basis elements is returned via (int *fdim).
@  'sym_opt' must be 0, 1 or -1.
@  sym_opt = 0:  the full lattice is calculated.
@  sym_opt = 1:  a basis for the symmetric matrices is calculated
@  sym_opt =-1:  a basis for the skewsymmetric matrices is calculated
@      
@---------------------------------------------------------------------------
@
\**************************************************************************/

matrix_TYP **
formspace (matrix_TYP **B, int Banz, int sym_opt, int *fdim)
{
  int i,j,k,l,m, no;
  int dim, dd, anz;
  matrix_TYP *X, *Xk, **E;
  int **XX, **Bi;
  int **pos, **sign;

  
  dim = B[0]->cols;
  if(B[0]->rows != dim)
  {
    printf("error in formspace: non-square matrix in group 'B'\n");
    exit(3);
  }
  for(i=1;i<Banz;i++)
  {
     if(B[i]->rows != dim || B[i]->cols != dim)
     {
       printf("error in formspace: different dimension of group elements\n");
       exit(3);
     }
  }
  if(sym_opt == 1)
    dd = (dim *(dim+1))/2;
  else if(sym_opt == -1)
    dd = (dim *(dim-1))/2;
  else if(sym_opt == 0)
    dd = dim * dim;
  else
  {
     printf("error in formspace: invalid sym_opt value %d\n", sym_opt);
     exit(2);
  }
  pos = (int **)xmalloc(dim *sizeof(int *));
  for(i=0;i<dim;i++)
  {
    pos[i] = (int *)xmalloc(dim *sizeof(int));
  }
  sign = (int **)xmalloc(dim *sizeof(int *));
  for(i=0;i<dim;i++)
  {
    sign[i] = (int *)xmalloc(dim *sizeof(int));
  }
  if(sym_opt == 0)
  {
     no = 0;
     for(i=0;i<dim;i++)
       for(j=0;j<dim;j++)
       { pos[i][j] = no; sign[i][j] = 1; no++;}
  }
  if(sym_opt == 1)
  {
     no = 0;
     for(i=0;i<dim;i++)
       for(j=0;j<=i;j++)
       { pos[i][j] = no; sign[i][j] = 1; no++;}
     no = 0;
     for(i=0;i<dim;i++)
       for(j=0;j<=i;j++)
       { pos[j][i] = no; sign[j][i] = 1; no++;}
  }
  if(sym_opt == -1)
  {
     no = 0;
     for(i=0;i<dim;i++)
       for(j=0;j<i;j++)
       { pos[i][j] = no; sign[i][j] = 1; no++;}
     no = 0;
     for(i=0;i<dim;i++)
       for(j=0;j<i;j++)
       { pos[j][i] = no; sign[j][i] = -1; no++;}
     for(i=0;i<dim;i++)
     { pos[i][i] = 0; sign[i][i] = 0;}
  }

  X = init_mat(dd * Banz, dd, "");
  XX = X->array.SZ;
  if(sym_opt == 0)
  {
    no = 0;
    for(i=0;i<Banz;i++)
    {
       Bi = B[i]->array.SZ;
       for(j=0;j<dim;j++)
         for(k=0;k<dim;k++)
         {
           for(l=0;l<dim;l++)
             for(m=0;m<dim;m++)
               XX[no][pos[l][m]] += sign[l][m] * Bi[l][j] * Bi[m][k];
           XX[no][pos[j][k]]--;
           no++;
         }
    }
  }
  if(sym_opt == 1)
  {
     no = 0;
     for(i=0;i<Banz;i++)
     {
       Bi = B[i]->array.SZ;
       for(j=0;j<dim;j++)
         for(k=0;k<=j;k++)
         {
           for(l=0;l<dim;l++)
             for(m=0;m<dim;m++)
               XX[no][pos[l][m]] += sign[l][m] * Bi[l][j] * Bi[m][k];
           XX[no][pos[j][k]]--;
           no++;
         }
     }
  }
  if(sym_opt == -1)
  {
     no = 0;
     for(i=0;i<Banz;i++)
     {
       Bi = B[i]->array.SZ;
       for(j=0;j<dim;j++)
         for(k=0;k<j;k++)
         {
           for(l=0;l<dim;l++)
             for(m=0;m<dim;m++)
               XX[no][pos[l][m]] += sign[l][m] * Bi[l][j] * Bi[m][k];
           XX[no][pos[j][k]]--;
           no++;
         }
     }
  }
  Xk = long_kernel_mat(X);
  free_mat(X);
  if(Xk == NULL)
   anz = 0;
  else
    anz = Xk->cols;
  if(anz != 0)
  {
    E = (matrix_TYP **)xmalloc(anz *sizeof(matrix_TYP *));
  }
  else
     E = NULL;
  for(i=0;i<anz;i++)
  {
    E[i] = init_mat(dim, dim, "");
    if(sym_opt == 0)
    {
       for(j=0;j<dim;j++)
         for(k=0;k<dim;k++)
            E[i]->array.SZ[j][k] = Xk->array.SZ[pos[j][k]][i];
    }
    if(sym_opt == 1)
    {
       for(j=0;j<dim;j++)
         for(k=0;k<=j;k++)
         {
            E[i]->array.SZ[j][k] = Xk->array.SZ[pos[j][k]][i];
            if(j != k)
              E[i]->array.SZ[k][j] = E[i]->array.SZ[j][k];
         }
    }
    if(sym_opt == -1)
    {
       for(j=0;j<dim;j++)
         for(k=0;k<j;k++)
         {
            E[i]->array.SZ[j][k] = Xk->array.SZ[pos[j][k]][i];
              E[i]->array.SZ[k][j] = -E[i]->array.SZ[j][k];
         }
       for(j=0;j<dim;j++)
        E[i]->array.SZ[j][j] = 0;
    }
  }
  for(i=0;i<anz;i++)
    Check_mat(E[i]);
  if(anz != 0)
    free_mat(Xk);
  for(i=0;i<dim;i++)
  { free(pos[i]); free(sign[i]);}
  free(pos); free(sign);
  *fdim = anz;
  for(i=0;i<anz;i++)
  {
    for(j=0;j<dim && E[i]->array.SZ[j][j] == 0; j++);
    if(j<dim && E[i]->array.SZ[j][j] < 0)
    {
       for(k=0;k<dim;k++)
        for(l=0;l<dim;l++)
          E[i]->array.SZ[k][l] = -E[i]->array.SZ[k][l];
    }
  }
  return(E);
}
