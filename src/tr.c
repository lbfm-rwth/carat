#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"

int main (int argc, char *argv[])
{
  int anz,
      i;

  char comment[1000];

  matrix_TYP **F,
              *tmp;

  read_header(argc, argv);

  if((FILEANZ != 1) || (is_option('h')))
  {
    printf("Usage: %s 'file' \n",argv[0]);
    printf("\n");
    printf("file: matrix_TYP.\n");
    printf("\n");
    printf(" Transposes the given matrices.\n");
    printf("\n");
    printf("cf. Tr_Bravais.\n");
    
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  F = mget_mat(FILENAMES[0],&anz);

  printf("#%d\n",anz);

  for (i=0;i<anz;i++){
     Check_mat(F[i]);
     tmp = tr_pose(F[i]);

     sprintf(comment,"tranposed of %d-th matrix of %s",i+1,FILENAMES[0]);

     put_mat(tmp,NULL,comment,2);
     free_mat(tmp);
  }

  exit(0);
}

