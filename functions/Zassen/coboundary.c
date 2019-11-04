#include "typedef.h"
#include "longtools.h"
#include "zass.h"
#include "matrix.h"


/**********************************************************************
@
@----------------------------------------------------------------------
@
@  static void coboundary(bravais_TYP *G,
@                        matrix_TYP *C,
@                        matrix_TYP *T)
@
@----------------------------------------------------------------------
@
***********************************************************************/

void 
coboundary (bravais_TYP *G, matrix_TYP *C, matrix_TYP *T)
{ 

  int i,
      k,
      l,
      n;

  matrix_TYP *A;
  matrix_TYP **erg;
  
  n = G->dim;

  A = init_mat(n*G->gen_no,n+1,"k");

  /* belegen der matrix A */
  for (i=0;i<G->gen_no;i++){
    for (k=0;k<n;k++){
      for (l=0;l<n;l++){
    	if (k == l){
          A->array.SZ[k+i*n][l] = G->gen[i]->array.SZ[k][l] - 1;
        }
        else{
          A->array.SZ[k+i*n][l] = G->gen[i]->array.SZ[k][l];
        }
      }
    }
  }
  
  for (i=0;i<G->gen_no;i++)
     for (k=0;k<n;k++)
        A->array.SZ[i*n+k][n] = C->array.SZ[i*n+k][0];

  k = long_row_gauss(A);
  real_mat(A,k,A->cols);

  /* kick out the rows which have 0's in the first n entries */
  for (i=A->rows-1;i>0;i--){
     for (k=0;k<n && A->array.SZ[i][k] == 0;k++);
     if (k == n){
        if (A->array.SZ[i][n] % C->kgv){
           fprintf(stderr,"error in coboundary\n");
           exit(3);
        }
        else{
           real_mat(A,A->rows-1,A->cols);
        }
     }
  }


  real_mat(C,A->rows,1);

  for (i=0;i<A->rows;i++)
     C->array.SZ[i][0] = A->array.SZ[i][n];

  real_mat(A,A->rows,n);

  Check_mat(A);
  Check_mat(C);

  /* put_mat(A,0,0,0);
  put_mat(C,0,0,0); */

  erg = long_solve_mat(A, C);

  if (erg == NULL || erg[0] == NULL){
     fprintf(stderr,"error in coboundary\n");
     exit(3);
  }

  /* put_mat(erg[0],0,0,0); */

  iscal_mul(T,erg[0]->kgv);
  for (i=0;i<G->dim;i++){
     T->array.SZ[i][G->dim] = erg[0]->array.SZ[i][0];
  }
  T->kgv = erg[0]->kgv;
  Check_mat(T);

  free_mat(A); A = NULL;
  for(i=0;i<2;i++){
     if (erg[i] != NULL) free_mat(erg[i]);
  }
  free(erg);

  return;
}


