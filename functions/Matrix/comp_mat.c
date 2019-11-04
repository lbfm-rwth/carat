#include "typedef.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: comp_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@
int cmp_mat(A, B)
matrix_TYP *A, *B;
@   compares the matrices A and B
@   if the first nonzero entry of A-B is positive (A>B) the result is 1
@   if the first nonzero entry of A-B is negative (A<B) the result is -1
@   if A == B the result is 0.
@
@   (the entries are assumed to be ordered in the following way:
@        (i,j)< (k,l) iff i<k or i==k and j<l
@---------------------------------------------------------------------------
@
\**************************************************************************/

int 
cmp_mat (matrix_TYP *A, matrix_TYP *B)
{
int i,j;
matrix_TYP *tmpA = NULL, *tmpB = NULL;

  if(A->cols != B->cols || A->rows != B->rows) {
     printf("matrices have incompatible dimensions\n");
     exit(3);
  }                    
  /* 
   *  Check_mat() is expensive, but our matrices are small
   */
  Check_mat( A );
  Check_mat( B );                    
  if ( A->array.N != NULL ) {
    if ( B->array.N == NULL ) {
      tmpB = copy_mat( B );
      kgv2rat( tmpB );
      tmpA = B; B = tmpB; tmpB = tmpA; tmpA = NULL;
    }
  } else {
    if ( B->array.N != NULL ) {
      tmpA = copy_mat( A );
      kgv2rat( tmpA );
      tmpB = A; A = tmpA; tmpA = tmpB; tmpB = NULL;
    }
  }
  if ( A->array.N == NULL ) {
    /*{{{}}}*/
    /*{{{  */
    if( A->kgv == B->kgv ) {
      for(i=0; i<A->rows; i++) {
        for(j=0; j<A->cols; j++) {
          if(A->array.SZ[i][j] != B->array.SZ[i][j]) {
            if(A->array.SZ[i][j] < B->array.SZ[i][j]) {
              return -1;
            } else {
              return  1;
            }
          }
        }
      }
    } else {
      for(i=0; i<A->rows; i++)
        for(j=0; j<A->cols; j++)
          if((A->array.SZ[i][j] * B->kgv) != (B->array.SZ[i][j] * A->kgv)) {
            if(A->array.SZ[i][j] < B->array.SZ[i][j]) {
              return -1;
            } else {
              return  1;
            }
          }
    }
    /*}}}  */
  } else {
    /*{{{  */
    for(i=0; i<A->rows; i++)
     for(j=0; j<A->cols; j++)
       if(    (A->array.SZ[i][j] != B->array.SZ[i][j] )
           || (A->array.N [i][j] != B->array.N [i][j] )
         ) {
           if ( A->array.SZ[i][j] == B->array.SZ[i][j] ) {
             if(A->array.N[i][j] > B->array.N[i][j]) {
               if ( tmpA ) {
                 free_mat( A );
               }
               if ( tmpB ) {
                 free_mat( B );
               }
               return -1;
             } else {
               if ( tmpA ) {
                 free_mat( A );
               }
               if ( tmpB ) {
                 free_mat( B );
               }
               return 1;
             }
           } else if( A->array.N[i][j] == B->array.N[i][j] ) {
             if(A->array.SZ[i][j] < B->array.SZ[i][j]) {
               if ( tmpA ) {
                 free_mat( A );
               }
               if ( tmpB ) {
                 free_mat( B );
               }
               return -1;
             } else {
               if ( tmpA ) {
                 free_mat( A );
               }
               if ( tmpB ) {
                 free_mat( B );
               }
               return  1;
             }
           } else { /* beware of overflow ! */
             if (   (A->array.SZ[i][j] * B->array.N[i][j])
                  < (B->array.SZ[i][j] * A->array.N[i][j])
                ) {
               if ( tmpA ) {
                 free_mat( A );
               }
               if ( tmpB ) {
                 free_mat( B );
               }
               return -1;
             } else {
               if ( tmpA ) {
                 free_mat( A );
               }
               if ( tmpB ) {
                 free_mat( B );
               }
               return  1;
             }
           }
         }
    /*}}}  */
  }
  return 0;
}
