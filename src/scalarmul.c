#include"typedef.h"
#include"tools.h"
#include"matrix.h"
#include"getput.h"
#include"longtools.h"



int main (int argc, char *argv[])
{
  matrix_TYP **M;

  int Manz,
         i,
         j,
         k,
         n,   /* the numerator */
         d,   /* the denominator */
         ggt;

  extern char **FILENAMES;

  char comment[1000];

  extern int FILEANZ;


  read_header(argc, argv);
  if((FILEANZ < 1) || (is_option('h'))){
     printf("Scalarmul file -n=numerator -d=denominator\n");
     printf("\n");
     printf("Multiplies all matrices in file with numerator/denominator.\n");
     printf("If either of these is ommitted, it is set to 1.\n");
     printf("\n");
     printf("\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }

  M = mget_mat(FILENAMES[0], &Manz);

  n = 1;
  if (is_option('n')){
     n = optionnumber('n');
  }

  d = 1;
  if (is_option('d')){
     d = optionnumber('d');
  }

  ggt = GGT(n,d);

  d = d/ggt;
  n = n/ggt;

  printf("#%d\n",Manz);

  for (i=0;i<Manz;i++){
     Check_mat(M[i]);
     rat2kgv(M[i]);
     if (n != 1){
        for (j=0;j<M[i]->rows;j++){
           for (k=0;k<M[i]->rows;k++){
              M[i]->array.SZ[j][k] = M[i]->array.SZ[j][k] * n;
           }
        }
     }
     M[i]->kgv = M[i]->kgv * d;
     Check_mat(M[i]);
     sprintf(comment,"%d-th matrix of %s multiplied by %d/%d",
                      i+1,FILENAMES[0],n,d);
     put_mat(M[i],NULL,comment,2);
     printf("\n");
  }

  exit(0);
}
