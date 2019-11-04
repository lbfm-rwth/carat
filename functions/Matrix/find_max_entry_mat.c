#include "typedef.h"
#include "matrix.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: find_max_entry.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

static int 
max_abs (int a, int b)
{ 
int aa = abs(a), ab=abs(b);

  return ( aa >= ab ? aa : ab );

}

/*{{{}}}*/
/*{{{  find_max_entry*/


/**************************************************************************\
@---------------------------------------------------------------------------
@ int find_max_entry(mat)
@ matrix_TYP *mat;
@
@ calculates the maximal entry of mat->array.SZ
@ if mat->array.Z != 0 or mat->prime != 0 the functions return 0.
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
find_max_entry (matrix_TYP *mat)
{
int i,j;
int **Z, *S_i;
boolean flag;

  if( mat->prime != 0 || mat->array.N != NULL ) {
    return 0;
  }
  
  flag = 0;
  Z = mat->array.SZ;
  if(mat->flags.Diagonal) {
    if(mat->flags.Scalar) {
      return abs(Z[0][0]);
    } else {
      for (i = 0; i < mat->rows; i++) {
        flag = max_abs(flag,Z[i][i]);
      }     
      return flag;
    }
  }
  if(mat->flags.Symmetric) {
    for(i = 0; i < mat->rows; i++) {
      for(j = i; j < mat->cols; j++) {
        if(Z[i][j]) flag = max_abs(flag,Z[i][j]);
      }
    } 
    return flag;
  }
  
  for(i = 0; i < mat->rows; i++) {
    S_i = Z[i];
    for(j = 0; j < mat->cols; j++) {
      if(S_i[j]) {
        if(abs(S_i[j]) > flag) {
          flag = abs(S_i[j]);
        }
      }
    }  
  }
  return flag;
}

/*}}}  */

