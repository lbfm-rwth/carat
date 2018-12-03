#include"typedef.h"
#include"tools.h"
#include"matrix.h"
#include"getput.h"


int main (int argc, char *argv[])
{
  matrix_TYP **M,
              *Trf;
  int erg,
      anz,
      i;

  read_header(argc, argv);
  if(FILEANZ != 1)
  {
     printf("Usage:  %s 'file' [-t]\n",argv[0]);
     printf("\n");
     printf("file:   matrix_TYP containing rational matrices.\n");
     printf("\n");
     printf("Calculates for each matrix in file a Gauss reduced matrix (row reduced).\n");
     printf("More precisely, by applying integral elementary transformations\n");
     printf("from the left, a staircase form of the input matrix is achieved.\n");
     printf("\n");
     printf("CAUTION: The program works with single precision and should be used only\n");
     printf("         for small examples.\n");
     printf("\n");
     printf("Options:\n");
     printf(" -t :     give the transforming matrices as well.\n");
     printf("\n");
     printf("Cf. Long_solve.\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }

  M = mget_mat(FILENAMES[0],&anz);
  for (i=0;i<anz;i++){
     Trf = init_mat(M[i]->rows, M[i]->rows, "");
     erg = Trf_gauss(M[i], Trf);
     put_mat(M[i], NULL, "row-gauss of matrix", 0);
     if (is_option('t')) put_mat(Trf, NULL, "Transformtion matrix", 0);
     free_mat(Trf);
     free_mat(M[i]);
  }

  exit(0);
}
