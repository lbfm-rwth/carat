#include "typedef.h"
#include "sort.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: compare.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
/*******************************************************************\
@ This file conatain functions that compare matrices or pointer
@ to integer with respect to different lexikographic orders.
@ The functions return -1 if the first arument is smaller with
@                         respect to the lexicographic order,
@                       0 if tboth arguments are equal
@                       1 if the first argument is bigger
@
@ If the arguments are 1-dimensional arrays, then the first
@ entry in which they differ decides.
@ For two dimensional arrays v < w if and only if
@    v[i][j] < w[i][j] and
@    v[k][j] == w[k][j] for all k < i and 0<=j<= dim
@    v[i][k] == w[i][k] for all k < j
@ For two dimensional symmetric arrays v < w if and only if
@    v[i][j] < w[i][j] and
@    v[k][j] == w[k][j] for all k < i and j<= i and
@    v[i][k] == w[i][k] for all k < j
@
@ for functions that compare first cols v < w if v^{tr} < w^{tr}
@ in the above sense.
\*******************************************************************/



/**************************************************************************\
@---------------------------------------------------------------------------
@ int mat_comp(m1, m2)
@ matrix_TYP *m1, *m2;
@
@ compares m1->array.SZ and m2->array.SZ
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
mat_comp (const matrix_TYP *m1, const matrix_TYP *m2)
{
  int i,
      j,
     *m1SZ,
     *m2SZ;

  if(m1->cols != m2->cols || m1->rows != m2->rows)
  {
    printf("cannot compare %dx%d with %dx%d-matrix\n", m1->rows, m1->cols, m2->rows, m2->cols);
    exit(3);
  }

  for(i=0;i<m1->rows;i++){
     m1SZ = m1->array.SZ[i];
     m2SZ = m2->array.SZ[i];
     for(j=0; j<m1->cols;j++,m1SZ++,m2SZ++){
       if(*m1SZ > *m2SZ)
          return(1);
       if(*m1SZ != *m2SZ)
          return(-1);
     }
  }

  return(0);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ int mat_col_comp(m1, m2)
@ matrix_TYP *m1, *m2;
@
@ compares m1->array.Sz and m2->array.SZ by going to the columns
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
mat_col_comp (const matrix_TYP *m1, const matrix_TYP *m2)
{
  int i,j;
  if(m1->cols != m2->cols || m1->rows != m2->rows)
  {
    printf("cannot compare %dx%d with %dx%d-matrix\n", m1->rows, m1->cols, m2->rows, m2->cols);
    exit(3);
  }
  for(i=0;i<m1->cols;i++)
  {
   for(j=0; j<m1->rows;j++)
   {
     if(m1->array.SZ[j][i] != m2->array.SZ[j][i])
     {
        if(m1->array.SZ[j][i] > m2->array.SZ[j][i])
           return(1);
        return(-1);
     }
   }
  }
   return(0);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ int lower_triangular_mat_comp(m1,m2)
@ matrix_TYP *m1, *m2;
@
@ compares the symmetric arrays m1->array.SZ and m2->array.SZ
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
lower_triangular_mat_comp (const matrix_TYP *m1, const matrix_TYP *m2)
{
  int i,j;
  if(m1->cols != m2->cols || m1->rows != m2->rows)
  {
    printf("cannot compare %dx%d with %dx%d-matrix\n", m1->rows, m1->cols, m2->rows, m2->cols);
    exit(3);
  }
  for(i=0;i<m1->rows;i++)
  {
   for(j=0; j<=i;j++)
   {
     if(m1->array.SZ[i][j] != m2->array.SZ[i][j])
     {
        if(m1->array.SZ[i][j] > m2->array.SZ[i][j])
           return(1);
        return(-1);
     }
   }
  }
   return(0);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ int vec_comp(v1, v2, dim)
@ int *v1, *v2, dim;
@
@ compares the 1-dimensional arrays v1 and v2 of length dim
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
vec_comp (const int *v1, const int *v2, int dim)
{
  int i;
  for(i=0;i<dim;i++)
  {
    if(v1[i] != v2[i])
    {
       if(v1[i] > v2[i])
         return(1);
       return(-1);
    }
  }
  return(0);
}




/**************************************************************************\
@---------------------------------------------------------------------------
@ int pointer_mat_comp(m1, m2, rows, cols)
@ int **m1, **m2, rows, cols;
@
@ compares the araays m1 and m2 of size rows x cols
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
pointer_mat_comp (int **m1, int **m2, int rows, int cols)
{
  int i,j;
  for(i=0;i<rows;i++)
  {
   for(j=0; j<cols;j++)
   {
     if(m1[i][j] != m2[i][j])
     {
        if(m1[i][j] > m2[i][j])
           return(1);
        return(-1);
     }
   }
  }
   return(0);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ int pointer_lower_triangular_mat_comp(m1,m2, n, m)
@ int **m1, **m2, n, m;
@
@ compares the symmetric arrays m1 and m2 of size n x m
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
pointer_lower_triangular_mat_comp (int **m1, int **m2, int n, int m)
{
  int i,j;
  for(i=0;i<n;i++)
  {
   for(j=0; j<=i;j++)
   {
     if(m1[i][j] != m2[i][j])
     {
        if(m1[i][j] > m2[i][j])
           return(1);
        return(-1);
     }
   }
  }
   return(0);
}
