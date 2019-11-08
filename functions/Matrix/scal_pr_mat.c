#include "typedef.h"
#include "matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: scal_pr_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/*{{{}}}*/
/*{{{  scal_pr*/
/*---------------------------------------------------------------*\
@-------------------------------------------------------------------------
@ matrix_TYP *scal_pr ( vectors, form, truth )
@ matrix_TYP *vectors, *form ;
@ boolean truth;
@
@ Calculates:      tr          
@         vectors  * form * vectors
@ vectors is supposed to be a matrix of row-vectors. 
@ form is supposed to be a symmetric bilinear form; if
@  form == NULL, the standard scalar product is used.
@
@ truth is an option to simplify the scalar-products for lattices 
@-------------------------------------------------------------------------
\*---------------------------------------------------------------*/
matrix_TYP *scal_pr (matrix_TYP *vectors, matrix_TYP *form, boolean truth )
{  
int  rV, cV, i, j, k;
int **V, **F, **S, *z, kgv = 0;
matrix_TYP *res;

    /*============================================================*\
    || Option to simplify the scalar-products for lattices        ||
    \*============================================================*/
  if(!truth) {
    tgauss(vectors);
  }
  
  V  = vectors->array.SZ;
  rV = vectors->rows;
  cV = vectors->cols;
  F  =   form->array.SZ;
  
  res = init_mat(rV,rV,"iks");
  if( save_null_mat(vectors) || save_null_mat(form)) {
    res->kgv = 0;
    res->flags.Integral  =
    res->flags.Symmetric =
    res->flags.Diagonal  = 
    res->flags.Scalar    = TRUE;
    return(res);
  }
  
  z = (int *)calloc(cV,sizeof(int));

  S = res->array.SZ;
  
    /*------------------------------------------------------------*\
    | Matrix product                         |
    | 5.02.92: because of the Gauss_Algorithm the vectors are      |
    | supposed to be given as a upper triangular matrix.           |
    \*------------------------------------------------------------*/
  for ( i = 0; i < rV; i++) {
      /*---------------------------------------------------------*\
      | Determine z as i-th row of vectors  * form        |
      \*---------------------------------------------------------*/
    for ( j = 0; j < cV; j++) {/* exchange '0' and 'i' */
      if ( V[i][j] != 0 ) {
        if (form->flags.Diagonal) {
          z[j] = V[i][j] * F[j][j];
        } else {
          for ( k = 0; k < cV; k++ ) {
            if ( F[j][k] != 0 ) {
              z[k] += V[i][j] * F[j][k];
            }                        
          }
        }                 
      }                        
    }
    /*---------------------------------------------------------*\
    |      tr              |
    | Determine i-th row of vectors  * form * vectors      |
    |      tr              |
    |  as z    * vectors          |
    \*---------------------------------------------------------*/
    for (j = 0 ; j < cV; j ++ ) {
      if ( z[j] != 0 )  {
        for ( k = 0; k <= i ; k++) {
          if ( V[k][j] != 0 ) {
            S[i][k] += z[j] * V[k][j]; 
          }
        }
        z[j] = 0;
      }
    }
  }
  
  for ( i = 0; i < rV; i++ ) {
    for ( j = 0; j < i; j++ ) {
      S[j][i] = S[i][j];      
    }
  }
  
  free(z);
  /*COUNTER--;*/
  kgv= vectors->kgv * vectors->kgv * form->kgv;
    /*------------------------------------------------------------*\
    | Fill data into matrix_TYP res                  |
    \*------------------------------------------------------------*/
  res->flags.Integral  = (kgv == 1);
  res->kgv     = kgv;
  Check_mat(res);
  return (res);
}
/*}}}  */
