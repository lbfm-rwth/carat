#include "typedef.h"
#include "matrix.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: gauss_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/*{{{}}}*/
/*{{{  tgauss*/
/*
@-------------------------------------------------------------------
@ int tgauss(mat)
@
@ Integral Gauss-algorithm where the matrix itself is changed
@ The Return-value is the rank of the matrix. The matrix is not
@ shrinked to its rank-size.
@-------------------------------------------------------------------
 */

int 
tgauss (matrix_TYP *mat)
{
int i,j,k,l;
int *nzero;
int **Z, *r, flag,ww;
int rg, rows, cols, tmp, actrow;

  /* put_mat(mat,NULL,NULL,0); */
  rows = mat->rows;
  cols = mat->cols;
  Z = mat->array.SZ;
  mat->flags.Symmetric = 
  mat->flags.Diagonal  = FALSE;
  nzero = (int *)malloc((cols+1)*sizeof(int));
  /*COUNTER++;*/
  actrow = 0;
  for (i = 0; i < cols; i++) {
    rg = 0;
    do {
      /* put_mat(mat,NULL,NULL,0); */
      /*
       *  Find the minimum of the non-zero entries of the i-th col
       */
      flag = 0;
      tmp = actrow;
      ww = ((actrow < rows) && (Z[actrow][i] != 0)) ? Z[actrow][i] : 0;
      for (j = actrow+1; j < rows; j++) {
        if (Z[j][i] != 0) {
          flag = 1;
          if(ww != 0) {
            if ( abs(Z[j][i])  < abs(ww)) {
              ww = Z[j][i];
              tmp = j;
            }
          } else {
            ww = Z[j][i];
            tmp = j;
          }
        }
      }
    /*===========================================================*\
    || Swap the tmp-row into the actrow-th rows and clear the    ||
    || i-th column downwards.                                    ||
    \*===========================================================*/
    if (ww != 0) rg = 1;
    if (flag) {
    if (tmp != actrow) {
      r = Z[actrow]; Z[actrow] = Z[tmp]; Z[tmp] = r;
      }
    k = 0;
    for(j = i; j < cols; j++) {
      if (Z[actrow][j] != 0) {
        nzero[k++] = j;
        }
      }
    nzero[k] = -1;
    for(tmp = actrow+1; tmp < rows; tmp++)
      if(Z[tmp][i] != 0) {
        k = 0;
        int waste= min_div(Z[tmp][i], Z[actrow][i]);
  /* printf("mindiv: %d\n",waste); */
        while(nzero[k] >= 0) {
          Z[tmp][nzero[k]] -= waste * Z[actrow][nzero[k]];
          k++;
          }
        }
      }
    /*
     * repeat the steps until every entry down from the i-th row 
     * is zero.
     */
    } while(flag);
    actrow += rg;
  }
  
  /*
   * determine the rank of the matrix                              
   */
  for(i = 0; i < cols; i++) {
    nzero[i] = 0;           
  }
  i = rows-1;
  tmp = 4*cols;
  while((i >= 0) && (memcmp(Z[i],nzero,tmp) == 0)) {
    i--;                                           
  }
  i++;
  rows = i;
  
    /*===========================================================*\
    || 03.02.92: Try to minimize the entries of the upper right  ||
    || triangular matrix by doing the subtraction upwards        ||
    \*===========================================================*/
  
  memset(nzero,0,cols*sizeof(int));
  for( i = 1; i < rows; i++) {
    tmp = 0;
    k = 0;
    while(Z[i][tmp] == 0) {
      tmp++;
    }
    for (j = tmp; j < cols; j++) {
      if(Z[i][j] != 0 ) {
        nzero[k++] = j; 
      }
    }
    for(j = i-1; j >= 0; j--) {
      int waste= min_div(Z[j][tmp],Z[i][tmp]);
      if (waste != 0) {
        for (l = 0; l <k; l++) {
          Z[j][nzero[l]] -= waste * Z[i][nzero[l]];
        }
   /*put_mat(mat,NULL,NULL,0);                    */
      }
    }
  }
  free(nzero);
  return(rows);
}





/**************************************************************************\
@---------------------------------------------------------------------------
@ int row_gauss(M)
@ matrix_TYP *M;
@
@ integral gauss applied to the rows of M,
$ return is the rank of the matrix.
@ M is changed!
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
row_gauss (matrix_TYP *M)
{
  int i,j;
  int step;
  int a1,a2,gcd;
  int tester, *v;
  int spos = 0;
  int x,y;

  for(step = 0; step < M->rows; step++)
  {
    tester = FALSE;
    while(tester == FALSE)
    {
       i = step;
       while(i<M->rows && M->array.SZ[i][spos] == 0)
          i++;
       if(i<M->rows)
       {
         tester = TRUE;
         if(i != step)
         {
           v = M->array.SZ[i];
           M->array.SZ[i] = M->array.SZ[step];
           M->array.SZ[step] = v;
         }
       }
       else
        spos++;
       if(spos == M->cols)
         return(step);
    }
    for(i=step+1;i<M->rows;i++)
    {
      if(M->array.SZ[i][spos] != 0)
      {
        gcd_darstell(M->array.SZ[step][spos], M->array.SZ[i][spos], &a1, &a2, &gcd);
       if(a1 == 0)
       {
         v = M->array.SZ[step];
         M->array.SZ[step] = M->array.SZ[i];
         M->array.SZ[i] = v;
         j = a1; a1 = a2; a2 = j;
       }
       M->array.SZ[step][spos] /= gcd;
       M->array.SZ[i][spos] /= gcd;
       for(j=spos + 1; j<M->cols; j++)
       {
         x = M->array.SZ[step][j];
         y = M->array.SZ[i][j];
         M->array.SZ[step][j] = a1 * x + a2 * y;
         M->array.SZ[i][j] = M->array.SZ[step][spos] * y - M->array.SZ[i][spos] * x;
       }
        M->array.SZ[step][spos] = gcd;
        M->array.SZ[i][spos] = 0;
      }
    }
  }
  return(M->rows);
}

/*}}}  */
/*{{{  ggauss*/
/*
@-------------------------------------------------------------------
@ matrix_TYP *ggauss(mat)
@ matrix_TYP *mat;
@
@ Integral Gauss-algorithm where a copy of mat is changed
@ and the number of rows is shrinked to the rank of the matrix.
@-------------------------------------------------------------------
 */
matrix_TYP *
ggauss (matrix_TYP *mat)
{
matrix_TYP *elt;
int rank;

  elt = copy_mat( mat );
/******************
  rank = tgauss( elt );
******************/
  rank = row_gauss( elt );
  if(rank != elt->rows) {
    real_mat( elt, rank, elt->cols);
  }
  return( elt );
}

/*}}}  */
