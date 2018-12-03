#include"typedef.h"
#include"getput.h"
#include"matrix.h"



int main (int argc, char *argv[])
{

	matrix_TYP **x, **S, *E;
        int anz1,
            anz2,
            i,
            j;

        char comment[1000];

        read_header(argc, argv);
        if(FILEANZ != 2)
        {
           printf("Usage: Mul 'file1' 'file2' [-x]\n");
           printf("\n");
           printf("file1: matrix_TYP containing the matrices A_i for 1<= i <= n.\n");
           printf("file2: matrix_TYP containing the matrices B_j for 1<= j <= m.\n");
           printf("\n");
           printf("Multiply the matrices in 'file1' with those in 'file2'.\n");
           printf("Called without any option the program ouptputs the products\n");
           printf("A_k * B_k for 1 <= k <= Min(n,m).\n");
           printf("\n");
           printf("Options:\n");
           printf(" -x    : perform all products of the form A_i * B_j for 1 <= i <=n,\n");
           printf("         1 <= j <= m.\n");
           printf("\n");
           printf("Cf. Add\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }
        x = mget_mat(FILENAMES[0], &anz1);
        S = mget_mat(FILENAMES[1], &anz2);

        if (is_option('x')){
           printf("#%d\n",anz1*anz2);
           for(i=0;i<anz1;i++){
              for(j=0;j<anz2;j++){
                 sprintf(comment,"product of the %d-th matrix of %s with the %d-th matrix of %s",i+1,FILENAMES[0],j+1,FILENAMES[1]);
                 E = mat_mul(x[i], S[j]);
                 put_mat(E,NULL,comment,2);
                 free_mat(E);
              }
           }
        }
        else{
           if (anz1>anz2) anz1 = anz2;
           printf("#%d\n",anz1);
           for(i=0;i<anz1;i++){
              sprintf(comment,"product of the %d-th matrix of %s with the %d-th matrix of %s",i+1,FILENAMES[0],i+1,FILENAMES[1]);
              E = mat_mul(x[i], S[i]);
              put_mat(E,NULL,comment,2);
              free_mat(E);
           }
        }

   exit(0);
}
