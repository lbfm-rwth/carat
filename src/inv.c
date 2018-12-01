#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"


int SFLAG;

int main (int argc, char *argv[])
{
  int anz,
      i;

  matrix_TYP **F,
              *tmp;

  read_header(argc, argv);
  if(FILEANZ != 1)
  {
    printf("usage: Inv file\n");
    printf(" where file contains a matrix_TYP.\n");
    printf("\n");
    printf(" Inverts the matrices given in file.\n");
    printf("\n");
    printf("\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  if (is_option('h') && optionnumber('h')==12){
     SFLAG = 1;
  }

  F = mget_mat(FILENAMES[0],&anz);

  printf("#%d\n",anz);

  for (i=0;i<anz;i++){
     Check_mat(F[i]);
     tmp = mat_inv(F[i]);
     rat2kgv(tmp);
     Check_mat(tmp);
     put_mat(tmp,NULL,NULL,2);
     free_mat(tmp);
     free_mat(F[i]);
  }

  free(F);
  if (is_option('h') && optionnumber('h')==12){
     pointer_statistics(0,0);
  }

  exit(0);
}

