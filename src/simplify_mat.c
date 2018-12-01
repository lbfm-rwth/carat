#include "typedef.h"
#include "getput.h"
#include "bravais.h"
#include "matrix.h"
#include "tools.h"

int gcd_mat(matrix_TYP *);

int main (int argc, char *argv[])
{
  int anz,
      i;

  rational *r;

  matrix_TYP **F;

  read_header(argc, argv);
  if(FILEANZ != 1)
  {
    printf("usage: Simplify_mat file\n");
    printf(" where file contains a matrix_TYP.\n");
    printf("\n");
    printf("Simplifies all the matrices in the given file, ie.\n");
    printf("divides all entiries of the matrix by the gcd of the entries.\n");
    printf("The gcd is taken to be positive!\n");
    printf("\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  F = mget_mat(FILENAMES[0],&anz);

  r = (rational *) malloc( sizeof(rational));

  printf("#%d\n",anz);

  for (i=0;i<anz;i++){
     Check_mat(F[i]);
     r->z = 1;
     r->n = abs(gcd_mat(F[i]));
     rscal_mul(F[i],*r); 
     put_mat(F[i],NULL,NULL,2);
  }

  exit(0);
}

int gcd_mat(matrix_TYP *A)
/* calculates the greatest common divisor of all entries in the
   matrix A */
{
  int i,
      j,
      erg;

  if (A->flags.Integral){
     erg = A->array.SZ[0][0];
     i=0;
     while ((i<A->rows) && (erg != 1) && (erg != (-1))){
        while ((j<A->rows) && (erg != 1) && (erg != (-1))){
           erg = GGT(erg,A->array.SZ[i][j]);
           j++;
        }
        i++;
     }
  }
  else{
     printf("matrix is not integral\n");
     exit(3);
  }

  return erg;
}
