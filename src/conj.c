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
             **Binv,
              *tmp,
              *tmp2;

  char comment[1000];

  read_header(argc, argv);
  if(FILEANZ < 2)
  {
     printf("Usage: %s 'file1' 'file2' [-x]\n",argv[0]);
     printf("\n");
     printf("file1: matrix_TYP containing matrices A_i\n");
     printf("file2: matrix_TYP containing matrices B_j\n");
     printf("\n");
     printf("Conjugates the i-th matrix of file1 with the i-th matrix  of\n");
     printf("file2, i.e. performs products of the form B_i * A_i * B_i^-1.\n");
     printf("\n");
     printf("Options:\n");
     printf("-x:    Conjugates all matrices of file1 with all matrices\n");
     printf("       of file2, i.e. performs all products of the for B_j * A_i * B_j^-1.\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }

  A = mget_mat(FILENAMES[0],&anz1);
  B = mget_mat(FILENAMES[1],&anz2);

  Binv = (matrix_TYP **) malloc(anz2 * sizeof(matrix_TYP *));

  for (i=0;i<anz2;i++){
     Binv[i] = mat_inv(B[i]);
  }

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
           tmp2 = mat_mul(A[i],Binv[j]);
           tmp = mat_mul(B[j],tmp2);
           free_mat(tmp2);
           sprintf(comment,
             "B*A*B^{-1} for A = %d-th matrix of %s, B= %d-th matrix of %s",
                            i+1,FILENAMES[0],j+1,FILENAMES[1]);
           put_mat(tmp,NULL,comment,2);
           free_mat(tmp);
        }
     }
  }
  else{
     printf("#%d\n",min);
     for (i=0;i<min;i++){
        tmp2 = mat_mul(A[i],Binv[i]);
        tmp = mat_mul(B[i],tmp2);
        free_mat(tmp2);
        sprintf(comment,
          "B*A*B^{-1} for A = %d-th matrix of %s, B= %d-th matrix of %s",
                         i+1,FILENAMES[0],i+1,FILENAMES[1]);
       put_mat(tmp,NULL,comment,2);
       free_mat(tmp);
     }
  }

  exit(0);
}
