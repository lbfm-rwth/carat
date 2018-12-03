#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"

int main (int argc, char *argv[])
{
  int anz,
      i;

  matrix_TYP **F;

  read_header(argc, argv);
  if(FILEANZ != 1)
  {
    printf("usage: Full file\n");
    printf(" where file contains a matrix_TYP.\n");
    printf("\n");
    printf(" Outputs the given matrices in full form, ie. without\n");
    printf(" abbreviation like 3d0, 3x0.\n");
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

  printf("#%d\n",anz);

  for (i=0;i<anz;i++){
     Check_mat(F[i]);
     put_mat(F[i],NULL,NULL,1);
  }

  exit(0);
}
