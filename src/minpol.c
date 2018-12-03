#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>

#include <longtools.h>
#include <bravais.h>
#include <idem.h>

int INFO_LEVEL;
extern int SFLAG;

int main(int argc,char **argv){

  matrix_TYP **x,
              *y,
              *z;

  char comment[1000];

  int i,
      j,
      anz;

  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) || (FILEANZ < 1)){
     printf("Usage: %s file [-r]\n",argv[0]);
     printf("\n");
     printf("file: matrix_TYP containing a set of integral square matrices.\n");
     printf("\n");
     printf("Calculates the minimal polynomial of each matrix in file.\n");
     printf("WARNING: If a given matrix is not integral, the result will be incorrect.\n");
     printf("\n");
     printf("-r    : calculate the integral roots of the given polynomials.\n");
     printf("\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }

  INFO_LEVEL = optionnumber('h');

  if (INFO_LEVEL & 12){
     SFLAG = 1;
  }

  x = mget_mat(FILENAMES[0],&anz);

  for (j=0;j<anz;j++){
     y = min_pol(x[j]);

     printf("The minimal polynomial to the %d matrix of %s\n",j+1,FILENAMES[0]);

     if (y->array.SZ[0][y->cols-1] == 1)
        printf(" X^%d",y->cols-1);
     else
        printf(" %d * X^%d",y->array.SZ[0][y->cols-1],y->cols-1);

     for (i=y->cols-2;i>0;i--)
        if (y->array.SZ[0][i] > 1)
           printf(" + %d * X^%d",y->array.SZ[0][i],i);
        else if (y->array.SZ[0][i] < -1)
           printf(" - %d * X^%d",-(y->array.SZ[0][i]),i);
        else if (y->array.SZ[0][i] == -1)
           printf(" - X^%d",i);
        else if (y->array.SZ[0][i] == 1)
           printf(" + X^%d",i);

     if (y->array.SZ[0][0]>0)
        printf(" + %d \n",y->array.SZ[0][0]);
     else if (y->array.SZ[0][0]<0)
        printf(" - %d \n",-(y->array.SZ[0][0]));

     printf("\n");
     /* put_mat(y,NULL,NULL,2); */

     if (is_option('r')){
        z = zeros(y);
        printf("Integral roots of the polynomial\n");
        put_mat(z,NULL,NULL,2);
        free_mat(z);
     }

     free_mat(y);
     free_mat(x[j]);
  }

  free(x);

  if (INFO_LEVEL & 12){
     pointer_statistics(0,0);
  }

  exit(1);

} /* main */
