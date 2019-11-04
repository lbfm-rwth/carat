#include "typedef.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: red_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

#if 0
/*{{{  */

static void 
sdec_mat (matrix_TYP *mat, matrix_TYP *trf)

{
int dim, i, j, flag;
int **M, **m1;

M = mat->array.SZ;
dim= mat->cols;
dim--;
flag = (M[dim][dim] == 0);

while(flag) {
  for(i = 0; flag && (i < dim); i++)
    flag = (M[dim][i] == 0);
  if (flag)     flag = ((--dim > 0) && (M[dim][dim] == 0));
  }
if (++dim < mat->cols) {
  for (i = 0; i < mat->cols; i++)
    if (i < dim)
      M[i] = (int *)realloc(M[i], dim*sizeof(int));
    else  {
      free(M[i]);
      }
  mat->array.SZ = (int **)realloc(M,dim*sizeof(int *));
  mat->cols = mat->rows = trf->rows = dim;
  }
m1 = (int **)malloc((dim+1) * sizeof(int *));
m1[dim] = (int *)malloc((dim) * sizeof(int ));
  for (i = 0; i < dim; i++)
           {
           m1[i]=mat->array.SZ[dim-1-i];
           }
  for (i = 0; i < dim; i++)
           {
           mat->array.SZ[i]=m1[i];
           }
  for (j = 0; j < dim; j++) {
  for (i = 0; i < dim; i++)
           m1[dim][i]=mat->array.SZ[j][dim-1-i];
  for (i = 0; i < dim; i++)
           mat->array.SZ[j][i]=m1[dim][i];
           }
free(m1[dim]);
free(m1);
}
/*}}}  */
#else
/*{{{   from ldec_mat*/
static void 
sdec_mat (matrix_TYP *mat, matrix_TYP *trf)
{
int j,dim, i, flag;
int **M;
int **m1, **t1;

  M = mat->array.SZ;
  dim = mat->cols;
  dim--;
  flag = (M[dim][dim] == 0);
  
  while (flag) {
    for (i = 0; flag && (i < dim); i++) {
      flag = (M[dim][i] == 0);
    }
    if (flag) {
      flag = ((--dim > 0) && (M[dim][dim] == 0));
    }
  }
  if (++dim < mat->cols) {
    for (i = dim; i < mat->cols; i++) {
      free(M[i]);
      free(trf->array.SZ[i]);
    }
    for (i = 0; i < dim; i++) {
      M[i] = (int *)realloc(M[i], dim* sizeof(int));
    }
    mat->array.SZ = (int **)realloc(M,dim*sizeof(int *));
    trf->array.SZ = (int **)realloc(trf->array.SZ,dim*sizeof(int *));
    trf->rows = mat->cols = mat->rows = dim;
  }
  m1 = (int **)malloc((dim+1) * sizeof(int *));
  m1[dim] = (int *)malloc((dim) * sizeof(int ));
  t1 = (int **)malloc(dim * sizeof(int *));
  for (i = 0; i < dim; i++) {
    m1[i]=mat->array.SZ[dim-1-i];
    t1[i]=trf->array.SZ[dim-1-i];
  }
  for (i = 0; i < dim; i++) {
    mat->array.SZ[i]=m1[i];
    trf->array.SZ[i]=t1[i];
  }
  for (j = 0; j < dim; j++) {
    for (i = 0; i < dim; i++) {
      m1[dim][i]=mat->array.SZ[j][dim-1-i];
    }
    for (i = 0; i < dim; i++) {
      mat->array.SZ[j][i]=m1[dim][i];
    }
  }
  free(m1[dim]);
  free(m1);
  free(t1);
}

/*}}}  */
#endif

/*
@
@ boolean extension(Mat,Trf);
@
@ Zusatzprogramm zu mat_red(): Ueberprueft bei M[i][i] = 2*M[i][j], ob
@ eine Addition zu weiterer Reduzierung fuehren kann.
@
static boolean extension(Mat,Trf)
matrix_TYP *Mat, *Trf;
{
int i,j,k, temp, ssignum;
int **M, *Z;
boolean flag = 0;

  M = Mat->array.SZ;
  Z = (int *)calloc(Mat->rows, sizeof(int));
  for(i = 1; i < Mat->rows; i++) {
    for(j = 0; (j < i) && (M[i][i] != 0); j++) {
      if(abs(2*M[i][j]) >= abs(M[i][i])) {
        if (M[i][i]*M[i][j] > 0) {
          ssignum = -1;
          for (k = 0;k < Mat->rows; k++) {
            Z[k] = M[j][k] - M[i][k];
          }
        } else {
          ssignum = 1;
          for(k = 0; k < Mat->rows; k++) {
            Z[k] = M[j][k] + M[i][k];
          }
        }
        flag = 1;
      }
      if(flag) {
        for (k = 0; k < Mat->rows; k++) {
          if ((k != i) && (k != j)) {
            if ( (temp =abs(2*Z[k])) > abs(M[k][k]) ||
                  temp               > abs(M[j][j]) ) {
              free(M[j]);
              M[j] = Z;
              if (ssignum > 0) {
                for (k = 0; k < Mat->rows; k++) {
                  M[k][j] += M[k][i];
                }
              } else {
                for (k = 0; k < Mat->rows; k++) {
                  M[k][j] -= M[k][i];
                }
              }
              return(TRUE);
            }
          }
        }
        flag = 0;
      }
    }
  }
  free(Z);
  return(FALSE);
}

*/

/*{{{}}}*/
/*{{{  smat_red*/
static matrix_TYP *
smat_red (matrix_TYP *mat)
{
int  **M, *M_i, *M_j, **m1;
rational *per, temp;
int i, j, k,l , flag;
int faktor, n;
int **T, *T_i, *T_j, **t1;
matrix_TYP *taM;

  n = mat->cols;
  taM = init_mat(n,n,"ik");
  M = mat->array.SZ;
  T = taM->array.SZ;
  
  per = (rational *)malloc(n*sizeof(rational));
  for ( i= 0; i< n; i ++) {
    T[i][i] =  1;
  }
  do {
    flag = 0;
    for(i = n-1; i >=0; i--) {
      M_i = M[i];
      T_i = T[i];
      for(j = n-1; j >=0; j--) {
        M_j = M[j];
        T_j = T[j];
        if (i != j) {
          if(M_i[i]) {
            if(abs(2*M_i[j]) > abs(M_i[i])) {
              faktor = (2*M_i[j])/M_i[i];
              faktor = (faktor / 2)+(faktor % 2);
              for(k = 0; k < n; k++) {
                M_j[k] -= faktor * M_i[k];
                T_j[k] -= faktor * T_i[k];
              }
              for(k = 0; k < n; k++) {
                M[k][j] -= faktor * M[k][i];
              }
              flag = 1;
            }
          } else {
            if (M_i[j]) {
              if ((faktor = M_j[j] / (M_i[j]*2)) != 0) {
                for(k = 0; k < n; k++) {
                  M_j[k] -= faktor * M_i[k];
                  T_j[k] -= faktor * T_i[k];
                }
                for (k = 0; k < n; k++) {
                  M[k][j] -= faktor*M[k][i];
                }
                flag = 1;
              }
              if(M_j[j] == 0) {
                for(k = 0; k < n; k++) {
                  if((i != k) && (M[k][j])) {
                    if((faktor = M[k][j] / M_i[j])!=0) {
                      for(l = 0; l < n; l++) {
                        M[k][j] -= faktor * M_i[l];
                        T[k][l] -= faktor * T_i[l];
                      }
                      for(l = 0; l < n ; l++) {
                        M[l][k] -= faktor * M[l][i];
                      }
                      flag = 1;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  } while ((flag) /* || extension(mat,taM) */ );
  /*
  * Sorting Mat in decreasing order of the absolute value of the
  * diagonal elements.
   */

  m1 = (int **)malloc((n+1) * sizeof(int *));
  t1 = (int **)malloc((n+1) * sizeof(int *));
  m1[n] = (int *)calloc(n , sizeof(int));
  
  for (i = n-1; i >= 0; i--) {
    per[i].z = abs(M[i][i]);
    per[i].n = i;
  }
  for( j = n-1; j >=1; j--) {
    temp = per[j];
    flag = j;
    for( i = j-1; i >= 0; i--) {
      if (per[i].z < temp.z) {
        temp = per[i];
        flag= i;
      }
    }
    if(flag != j) {
      per[flag] = per[j];
      per[j] = temp;
    }
  }
  i = n;
  while(per[--i].z == 0);
  ++i;
  while((i<n)&&(flag = memcmp(m1[n],M[per[i].n],n*sizeof(int)))) {
    i++;
  }
  if (i < n-1) {
    for(j = i+1; j < n; j++) {
      if( (flag = memcmp(m1[n],M[per[j].n],n*sizeof(int)))  ) {
        temp = per[i]; per[i] = per[j]; per[j] = temp;
        i++; j = i+1;
      }
    }
  }
  for (i = 0; i < n; i++) {
    m1[i] = M[i];
    t1[i] = T[i];
  }
  for(i = 0; i < n; i++) {
    M[i] = m1[per[i].n];
    T[i] = t1[per[i].n];
    memcpy(m1[n],M[i],n*sizeof(int));
    for(j = 0; j < n; j++) {
      M[i][j] = m1[n][per[j].n];
    }
  }
  free(m1[n]);
  free(m1);
  free(t1);
  free(per);
  
  return (taM);
}
/*}}}  */
/*{{{  sdec_mat*/

/*}}}  */
/*{{{  mat_red*/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *mat_red(Mat)
@ matrix_TYP *Mat;
@
@ applies a pair_reduction to the symmetric matrix Mat
@ and returns the transformation matrix
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
mat_red (matrix_TYP *Mat)
{
  return smat_red(Mat);
}

void 
dec_mat (matrix_TYP *Mat, matrix_TYP *Trf)
{
  sdec_mat(Mat,Trf);
}
/*}}}  */
