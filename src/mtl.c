#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"

int main (int argc, char *argv[])
{
  int anz,
      i,
      j,
      k,
      rows,
      cols;

  char comment[1000];

  matrix_TYP **F,
              *tmp;

  read_header(argc, argv);
  if(FILEANZ < 1 || is_option('h'))
  {
    printf("Usage: Mtl file\n");
    printf("\n");
    printf(" where file contains a matrix_TYP.\n");
    printf("\n");
    printf("Writes the entries of the matrices given in file in ONE\n");
    printf("matrix, which rows consists of the concatenation of the\n");
    printf("rows of the individual matrices in file.\n");
    printf("Note that this function does exactly the inverse of Ltm!\n"); 
    printf("\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  F = mget_mat(FILENAMES[0],&anz);

  rows = F[0]->rows;
  cols = F[0]->cols;

  tmp = init_mat(anz,rows*cols,"r");

  for (i=0;i<anz;i++){
     rat2kgv(F[i]);
     if (rows != F[i]->rows){
       printf("The given matrices must all have the same number of rows.\n");
       exit(3);
     }
     if (cols != F[i]->cols){
       printf("The given matrices must all have the same number of columns.\n");
       exit(3);
     }
     for (j=0;j<cols;j++){
        for (k=0;k<rows;k++){
           tmp->array.SZ[i][cols*k+j] = F[i]->array.SZ[k][j];
           tmp->array.N[i][cols*k+j] = F[i]->kgv;
        }
     }
  }

  rat2kgv(tmp);
  Check_mat(tmp);

  sprintf(comment,"Mtl on file %s",FILENAMES[0]);

  put_mat(tmp,NULL,comment,2);


  exit(0);
}

