#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"

main (int argc, char *argv[])
{
  int i,
      j,
      k,
      anz;

  matrix_TYP **A;


  read_header(argc, argv);
  if(FILEANZ < 1)
  {
    printf("Usage: Conv 'file' -g -G -m -M\n");
    printf("\n");
    printf("file: matrix_TYP.\n");
    printf("\n");
    printf("Converts the matrix_TYP into the format of other\n");
    printf("programs available. The options specify the program,\n");
    printf("and only one of them is allowedat one time.\n");
    printf("\n");
    printf("-g or -G: GAP-format\n");
    printf("-m or -M: MAPLE-format\n");
    printf("\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  A = mget_mat(FILENAMES[0],&anz);

  if (is_option('m') || is_option('M')){
    printf("with(linalg);\n");
    for (i=0;i<anz;i++){
       if (!A[i]->flags.Integral){
          rat2kgv(A[i]);
       }
       printf("a%d := array(1..%d,1..%d,[\n",i+1,A[i]->rows,A[i]->cols);
       for (j=0;j<A[i]->rows;j++){
          printf("[");
          for (k=0;k<A[i]->cols;k++){
             printf(" %d ",A[i]->array.SZ[j][k]);
             if (!A[i]->flags.Integral){
                printf("/%d",A[i]->kgv);
             }
             if (k!=(A[i]->cols-1)){
                printf(",");
             }
          }
          if (j != (A[i]->rows-1)){
             printf("],\n");
          }
          else{
             printf("]\n");
          }
       }
       printf("]);\n");
    }
  }
  else if(is_option('g') || is_option('G')){
    for (i=0;i<anz;i++){
       if (!A[i]->flags.Integral){
          rat2kgv(A[i]);
       }
       printf("a%d := [",i+1);
       for (j=0;j<A[i]->rows;j++){
          printf("[");
          for (k=0;k<A[i]->cols;k++){
             printf(" %d ",A[i]->array.SZ[j][k]);
             if (!A[i]->flags.Integral){
                printf("/%d",A[i]->kgv);
             }
             if (k!=(A[i]->cols-1)){
                printf(",");
             }
          }
          if (j != (A[i]->rows-1)){
             printf("],\n");
          }
          else{
             printf("]");
          }
       }
       printf("];\n");
    }
  }
  else{
    printf("you must specify at least one option. try -h first\n");
  }

  exit(0);
}
