#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"

int main (int argc, char *argv[])
{
  int anz,
      anz2,
      i;

  char comment[1000];

  matrix_TYP **Formen,
             **A,
              *tmp;

  read_header(argc, argv);
  if(FILEANZ < 2 || is_option('h'))
  {
    printf("Usage: Normlin file1 file2\n");
    printf("\n");
    printf("where file1 and file2 contain a matrix_TYP.\n");
    printf("\n");
    printf("Calculates for each matrix A in 'file2' a matrix X with the\n");
    printf("property that\n");
    printf(" \\sum_j X_{i,j} F_j = A^{tr} F_j A with F_j in 'file1'\n");
    printf("\n");
    printf("CAUTION: The matrix describes the action on rows!\n");
    printf("\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  Formen = mget_mat(FILENAMES[0],&anz);
  A = mget_mat(FILENAMES[1],&anz2);

  for (i=0;i<anz2;i++){
     tmp = normlin(Formen,A[i],anz);
     sprintf(comment,"Normalin for %d-th matrix of %s on the formspace %s",
                      i+1,FILENAMES[1],FILENAMES[2]);
     Check_mat(tmp);
     put_mat(tmp,NULL,comment,2);
     free_mat(tmp);
  }

  exit(0);
}
