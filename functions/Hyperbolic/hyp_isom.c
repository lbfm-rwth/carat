#include "typedef.h"
#include "matrix.h"

/************************************************************************** \
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: hyp_isom.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/



/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *hyperbolic_isometry(x1, x2, S)
@ matrix_TYP *x1, *x2, *S;
@
@ The first row in the matrix x1 (resp. x2) is interpreted as a vector
@ v1 (resp. v2) in Z^n.
@ The function calculates a matrix A in GL_n(Q) (if exists) with
@          v1A = v2 and ASA^{tr} = S
@ Any integral matrix with this property must be of the form BA where
@ B is a matrix of the group that can be caluclated with
@    B = hyperbolic_stabilizer(x1, S)
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *hyperbolic_isometry(x1, x2, S)
matrix_TYP *x1, *x2, *S;
{
  matrix_TYP *x1S, *x2S, *x1SO, *x2SO, **S1red, **S2red, *T1,*T1i, *T2;
  matrix_TYP *Y1, *Y1i, *Y2, *Y;
  int i,j,k, n;

  extern matrix_TYP *solve_mat();
  extern matrix_TYP *mat_mul();
  extern matrix_TYP *mat_inv();
  extern matrix_TYP *scal_pr();
  extern matrix_TYP *pr_isom();
  extern matrix_TYP *init_mat();

  n = x1->cols;
  if((S1red = (matrix_TYP **)malloc(1 *sizeof(matrix_TYP *))) == NULL)
  {
    printf("malloc of 'S1red' in 'hyperbolic_isometry' failed\n");
    exit(2);
  }
  if((S2red = (matrix_TYP **)malloc(1 *sizeof(matrix_TYP *))) == NULL)
  {
    printf("malloc of 'S2red' in 'hyperbolic_isometry' failed\n");
    exit(2);
  }
  x1S = mat_mul(x1, S);
  x1SO = solve_mat(x1S);
  S1red[0] = scal_pr(x1SO, S, 1);
  x2S = mat_mul(x2, S);
  x2SO = solve_mat(x2S);
  S2red[0] = scal_pr(x2SO, S, 1);
  Y1 = pr_isom(S1red, S2red, 1, NULL, 0, NULL);
  free_mat(S1red[0]); free(S1red);
  free_mat(S2red[0]); free(S2red);
  free_mat(x1S); free_mat(x2S);
  if(Y1 == NULL)
  {
    free_mat(x1SO); free_mat(x2SO);
    return(NULL);
  }
  T1 = init_mat(n,n,"");
  for(i=0;i<n-1;i++)
   for(j=0;j<n;j++)
     T1->array.SZ[i][j] = x1SO->array.SZ[i][j];
  for(i=0;i<n;i++)
   T1->array.SZ[n-1][i] = x1->array.SZ[0][i];
  T1i = mat_inv(T1);
  free_mat(T1);
  free_mat(x1SO);

  T2 = init_mat(n,n,"");
  Y1i = mat_inv(Y1);
  Y2 = mat_mul(Y1i, x2SO);
  for(i=0;i<n-1;i++)
    for(j=0;j<n;j++)
      T2->array.SZ[i][j] = Y2->array.SZ[i][j];
  for(i=0;i<n;i++)
    T2->array.SZ[n-1][i] = x2->array.SZ[0][i];
  free_mat(x2SO);
  free_mat(Y2);
  free_mat(Y1);
  free_mat(Y1i);

  Y = mat_mul(T1i, T2);
  free_mat(T1i); free_mat(T2);
  return(Y);
}
