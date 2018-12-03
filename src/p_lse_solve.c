#include "typedef.h"
#include "tools.h"
#include "matrix.h"
#include "getput.h"

int main (int argc, char *argv[])
{
  matrix_TYP **A, **B, **X;
  int i, prime, anz;
  int Aanz, Banz;

  extern char **FILENAMES;
  extern int FILEANZ;

  extern matrix_TYP **p_lse_solve();
  extern matrix_TYP **mget_mat();

  read_header(argc, argv);
  if(FILEANZ <= 0)
  {
     printf("usage: P_lse_solve 'file1' ['file2'] -p=prime\n");
     printf("file1: matrix_TYP with matrix A\n");
     printf("file2: matrix_TYP with matrix B\n");
     printf("\n");
     printf("If 'file2' is not given, B = 0\n");
     printf("Solves the linear equations\n");
     printf(" A * X = B\n");
     printf("modulo the given prime (default 1949).\n");
     printf("\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }
  A = mget_mat(FILENAMES[0], &Aanz);
  if(FILEANZ > 1)
    B = mget_mat(FILENAMES[1], &Banz);
  else
  {
    B = (matrix_TYP **) malloc(1 *sizeof(matrix_TYP *));
    B[0] = NULL;
  }
  prime = optionnumber('p');
  if(prime == 0)
   prime = 101;
  
  X = p_lse_solve(A[0], B[0], &anz, prime);
  if(B[0] != NULL){
     put_mat(X[0], NULL, "solution of the inhomogenous equation", 0);
     free_mat(X[0]);
  }

  for(i=1;i<anz;i++){
     put_mat(X[i], NULL, "solution of the homogenous equation", 0); 
     free_mat(X[i]);
  }
  if (X!=NULL) free(X);

  exit(0);
}
