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
     printf("Usage: %s 'file' [-c=n1] [-r=n2]\n",argv[0]);
     printf("\n");
     printf("file: matrix_TYP.\n");
     printf("\n");
     printf("Writes the lines of the matrices given in file in\n");
     printf("seperate matrices. The default assumption is that\n");
     printf("the resulting matrices should be square.\n");
     printf("In case either option r or c are given, the matrices\n");
     printf("are given as n2xn1 matrices, where n1*n2=#columns of the\n");
     printf("matrices of file.\n");
     printf("Note that this function does exactly the inverse of Mtl!\n");
     printf("\n");
     printf("OPTIONS:\n");
     printf("\n");
     printf("-c=n1:  the matrices will have n1 columns.\n");
     printf("-r=n2:  the matrices will have n2 rows.\n");
     printf("\n");
     printf("Cf. Mtl\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }

  F = mget_mat(FILENAMES[0],&anz);

  for (i=0;i<anz;i++){
     rat2kgv(F[i]);
     if ((F[i]->cols != F[0]->cols) || (F[i]->rows != F[0]->rows)){
        printf("The matrices given in %s must have same dimensions\n",
                 FILENAMES[0]);
        exit(3);
     }
  }

  if (is_option('r') && is_option('c')){
     rows=optionnumber('r');
     cols=optionnumber('c');
  }
  else if(is_option('r')){
     rows = optionnumber('r');
     cols = floor(F[0]->cols/rows + 0.5);
  }
  else if(is_option('c')){
     cols = optionnumber('c');
     rows = floor(F[0]->cols/cols + 0.5);
  }
  else{
     cols = floor(sqrt(F[0]->cols)+0.5);
     rows = floor(sqrt(F[0]->cols)+0.5);
  }

  if ((rows*cols) !=  F[0]->cols){
     printf("The specification is not suitable\n");
     exit(3);
  }

  printf("#%d\n",anz*F[0]->rows);

  tmp = init_mat(rows,cols,"");

  for (i=0;i<anz;i++){
     for (j=0;j<F[0]->rows;j++){
        tmp->kgv = F[i]->kgv;
        for  (k=0;k<F[0]->cols;k++){
           tmp->array.SZ[(k - (k%cols))/cols][k % cols]=F[i]->array.SZ[j][k];
        }
        sprintf(comment,"Ltm on file %s",FILENAMES[0]);
        Check_mat(tmp);
        put_mat(tmp,NULL,comment,2);
     }
  }


  exit(0);
}

