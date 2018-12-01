#include "typedef.h"
#include "getput.h"
#include "bravais.h"
#include "matrix.h"
#include "datei.h"
#include "tools.h"

int main (int argc, char *argv[])
{
  int anz,
      i,
      j,
      ggt,
      sum;

  char comment[1000];

  matrix_TYP **F;

  read_header(argc, argv);

  if((FILEANZ != 1) || (is_option('h')))
  {
    printf("Usage: %s 'file'\n",argv[0]);
    printf("\n");
    printf("file: matrix_TYP\n");
    printf("\n");
    printf("Calculates the trace of the given matrices.\n");
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
     rat2kgv(F[i]);

     if (F[i]->rows != F[i]->cols){
        printf("Trace of a non square matrix?\n");
        exit(3);
     }
     else{
       sum = 0;
       for (j=0;j<F[i]->cols;j++){
         sum = sum + F[i]->array.SZ[j][j];
       }
       if (!F[i]->flags.Integral){
         ggt = GGT(sum,F[i]->kgv);
         if (abs(F[i]->kgv)==ggt){
            printf("Trace of the %d-th matrix in %s: %d\n",
                i+1,FILENAMES[0],sum/ggt*(F[i]->kgv/ggt));
         }
         else{
            printf("Trace of the %d-th matrix in %s: %d/%d\n",
                 i+1,FILENAMES[0],sum/ggt,F[i]->kgv/ggt);
         }
       }
       else{
         printf("Trace of the %d-th matrix in %s: %d\n",i+1,FILENAMES[0],sum);
       }
     }
  }

  exit(0);
}
