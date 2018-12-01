#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"

int main (int argc, char *argv[])
{
  int i,
      j,
      anz1,
      anz2,
      min;

  matrix_TYP **A,
             **B,
              *tmp;

  char comment[1000];

  read_header(argc, argv);
  if((FILEANZ < 2) || (is_option('h')))
  {
    printf("Usage: Kron 'file1' 'file2' [-x]\n");
    printf("\n");
    printf("file1: matrix_TYP.\n");
    printf("file2: matrix_TYP.\n");
    printf("\n");
    printf("Calculates the Kronecker product of the matrices given\n");
    printf("in file1 with those of file2.\n");
    printf("\n");
    printf("OPTIONS:\n");
    printf("-x   : Outputs every possible Kronecker product. If this option\n");
    printf("       is not present, the Kronecker product of the i-th matrix\n");
    printf("       of file1 with the i-th matrix of file2 is calculated.\n");
    printf("\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  A = mget_mat(FILENAMES[0],&anz1);
  B = mget_mat(FILENAMES[1],&anz2);

  if (anz1<anz2){
     min = anz1;
  }
  else{
     min = anz2;
  }

  if (is_option('x')){
     printf("#%d\n",anz1*anz2);
     for (i=0;i<anz1;i++){
        for (j=0;j<anz2;j++){
           tmp = kron_mat(A[j],B[i]);
           sprintf(comment,
              "Kronecker product of %d-th matrix of %s with %d-th matrix of %s",
               j+1,FILENAMES[0],i+1,FILENAMES[1]);
           put_mat(tmp,NULL,comment,2);
           free_mat(tmp);
        }
     }
  }
  else{
     printf("#%d\n",min);
     for (i=0;i<min;i++){
        tmp = kron_mat(A[i],B[i]);
        sprintf(comment,
           "Kronecker product of %d-th matrix of %s with %d-th matrix of %s",
            i+1,FILENAMES[0],i+1,FILENAMES[1]);
        put_mat(tmp,NULL,comment,2);
        free_mat(tmp);
     }
  }


  exit(0);
}
