#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"

int main (argc,argv)
int argc;
char *argv[];
{
  int i,
      j,
      anz1,
      anz2,
      min,
      l,
      r,
      z,
      n;

  matrix_TYP **A,
             **B,
              *tmp;

  char comment[1000];

  rational l_rat,
           r_rat;


  read_header(argc, argv);
  if(FILEANZ < 2)
  {
    printf("usage: Add 'file1' 'file2' [-x] -l=n1 -r=n2 -n=n3\n");
    printf(" where file1 and file2 contain a matrix_TYP .\n");
    printf("\n");
    printf(" Calculates the sums (n1 * A + n2 * B)/n3 with matrices A, B\n");
    printf(" taken from file1, file2 respectively.\n");
    printf("\n");
    printf("-x:    Calculates all possible sums of matrices of file1\n");
    printf("       with file2.\n");
    printf("       If this option is not present, the programm adds the\n");
    printf("       i-th matrix of file1 to the i-th matrix of file2\n");
    printf("\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  A = mget_mat(FILENAMES[0],&anz1);
  B = mget_mat(FILENAMES[1],&anz2);

  if (anz1<anz2){
     min = anz1;
  }
  else{
     min = anz2;
  }

  if (is_option('n')){
     n = optionnumber('n');
  }
  else{
     n = 1;
  }

  if (is_option('l')){
     l = optionnumber('l');
  }
  else{
     l = 1;
  }

  if (is_option('r')){
     r = optionnumber('r');
  }
  else{
     r = 1;
  }

  l_rat.z = l;
  l_rat.n = n;

  r_rat.z = r;
  r_rat.n = n;

  if (is_option('x')){
     printf("#%d\n",anz1*anz2);
     for (i=0;i<anz1;i++){
        for (j=0;j<anz2;j++){
           tmp = mat_add(A[i],B[j],l_rat,r_rat);
           sprintf(comment,"sum of %d-th matrix of %s with %d-th matrix of %s",
                            i+1,FILENAMES[0],j+1,FILENAMES[1]);
           put_mat(tmp,NULL,comment,2);
           free_mat(tmp);
        }
     }
  }
  else{
     printf("#%d\n",min);
     for (i=0;i<min;i++){
       tmp = mat_add(A[i],B[i],l_rat,r_rat);
       sprintf(comment,"sum of %d-th matrix of %s with %d-th matrix of %s",
                        i+1,FILENAMES[0],i+1,FILENAMES[1]);
       put_mat(tmp,NULL,comment,2);
       free_mat(tmp);
     }
  }
  return 0;
}
