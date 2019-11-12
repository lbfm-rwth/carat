#include "typedef.h"
#include "utils.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: p_mat_det.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ int p_mat_det(Mat, prime)
@ matrix_TYP *Mat;
@ int prime;
@
@ Calculates the determinant of Mat modulo prime.
@ The entries of Mat are not changed.
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
p_mat_det (matrix_TYP *Mat, int prime)
{
  int i,j;
  int step, dim;
  int **M, *tmp;
  int det, msi;

  dim = Mat->cols;
  if(prime < 0)
    prime = -prime;
  if(Mat->rows != Mat->cols)
  {
    fprintf(stderr, "error: can't calculate the determinante of a non-square matrix\n");
    exit(3);
  }
  M = (int **)xmalloc(dim *sizeof(int *));
  for(i=0;i<dim;i++)
  {
      M[i] = (int *)xmalloc(dim *sizeof(int));
  }
  for(i=0;i<dim;i++)
    for(j=0;j<dim;j++)
      M[i][j] = Mat->array.SZ[i][j] % prime;
  det = 1;
  for(step = 0; step < dim;step++)
  {
    /************************************************************************\
    | search for non zero entry in the colmn no. step and
    | swap it to M[step][step]
    \************************************************************************/
     for(i=step; i<dim && M[i][step] == 0; i++);
     if(i == dim)
     {
         for(j=0;j<dim;j++)
          free(M[j]);
         free(M);
         return(0);
     }
     if(i != step){
        tmp = M[step]; M[step] = M[i]; M[i] = tmp; det = -det;
     }
     msi = p_inv(M[step][step], prime);
     i = (msi * M[step][step])%prime;
     if(i<0) i+=prime;
     if(i != 1)
     printf("TEST %d\n", i);

     /* multiplying the step-th row with msi */
     for(i=step+1;i<dim;i++){
        M[step][i] *= msi;
        M[step][i] %=prime;
     }

     /* changed tilman 9/12/96 from: 
     det *= msi; det %=prime;
     to :*/

     det *= M[step][step]; det %=prime;

     M[step][step] = 1;

    /************************************************************************\
    | clear column no. step
    \************************************************************************/
     for(i=step+1;i<dim;i++)
     {
       if(M[i][step] != 0)
       {
          for(j=step+1;j<dim;j++)
          {
            M[i][j] -= M[i][step] * M[step][j];
            M[i][j] %= prime;
          }
          M[i][step] = 0;
       }
     }
  }
  for(i=0;i<dim;i++)
     free(M[i]);
  free(M);
  det %= prime;
  if(det < -(prime/2))
     det += prime;
  if(det > (prime/2))
     det -=prime;
  if(prime == 2 && det == -1)
     det = 1;
  return(det);
}
