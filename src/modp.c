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
      prime;

  matrix_TYP **F,
              *tmp;

  read_header(argc, argv);
  if((FILEANZ != 1) || (!is_option('p')))
  {
    printf("usage: Modp file -p=prime\n");
    printf(" where file contains a matrix_TYP.\n");
    printf("\n");
    printf(" Writes the given matrices in file modulo the given prime.\n");
    printf("\n");
    printf("\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  F = mget_mat(FILENAMES[0],&anz);

  prime = optionnumber('p');

  printf("#%d\n",anz);

  for (i=0;i<anz;i++){
     Check_mat(F[i]);
     if (!F[i]->flags.Integral){
       printf("Modp uses only integral matrices \n");
       exit(3);
     }
     for (j=0;j<F[i]->rows;j++){
       for (k=0;k<F[i]->cols;k++){
       F[i]->array.SZ[j][k] = F[i]->array.SZ[j][k] % prime;
       }
     }
     Check_mat(F[i]);
     put_mat(F[i],NULL,NULL,2);
  }

  exit(0);
}

