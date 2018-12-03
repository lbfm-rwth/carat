#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"

int main (int argc, char *argv[])
{
  int anz,
      i,
      prime = 1949,
      det;

  matrix_TYP **F,
              *tmp;

  read_header(argc, argv);
  if(FILEANZ < 1 || is_option('h'))
  {
    printf("usage: Pdet file [-p=prime]\n");
    printf(" where file contains a matrix_TYP.\n");
    printf("\n");
    printf("Calculates the determinant of the matrices given in file\n");
    printf("modulo a prime. The default prime is 1949.\n");
    printf("\n");
    printf("Options:\n");
    printf("\n");
    printf("-p=prime : the algorithm is done modulo prime.\n");
    printf("\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  F = mget_mat(FILENAMES[0],&anz);

  if (is_option('p')){
     prime = optionnumber('p');
  }

  for (i=0;i<anz;i++){
     Check_mat(F[i]);
     det = p_mat_det(F[i],prime);
     printf("Determinant of the %d-th matrix in %s (modulo %d): %d \n",
             i+1,FILENAMES[0],prime,det);

  }

  exit(0);
}
