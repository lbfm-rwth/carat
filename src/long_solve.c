#include "typedef.h"
#include "getput.h"
#include "longtools.h"
#include "matrix.h"

int INFO_LEVEL;
int SFLAG;


int main (int argc, char *argv[])
{

	matrix_TYP **A, **B, **X;
        int i, Aanz, Banz;


        read_header(argc, argv);
        if(FILEANZ != 1 && FILEANZ != 2)
        {
           printf("Usage: %s 'file1' ['file2']\n", argv[0]);
           printf("\n");
           printf("file1: matrix_TYP containing the matrices A_i\n");
           printf("file2: matrix_TYP containing the matrices B_i\n");
           printf("\n");
           printf("Calculates matrices X_i such that A_i X_i = B_i for all matrices\n");
           printf("A_i in file1. If 'file2' is omitted the B_i are assumed to be 0.\n");
           printf("\n");
           printf("Cf. Gauss.\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }

        if (is_option('h')){ INFO_LEVEL = optionnumber('h');}
        if (INFO_LEVEL & 8){
           SFLAG = 1;
        }

	A = mget_mat (FILENAMES[0], &Aanz);
        if(FILEANZ == 2){
           B = mget_mat(FILENAMES[1], &Banz);
        }
        else{
           Banz = 0;
        }
        if(Aanz > 1)
          printf("#%d\n", Aanz);
        for(i=0;i<Aanz;i++)
        {
            if(i<Banz)
              X = long_solve_mat(A[i], B[i]);
            else
              X = long_solve_mat(A[i], NULL);
            if(X[0] != 0)
            {
               put_mat(X[0], NULL, "inhomogenous solution", 2);
               free_mat(X[0]);
            }
            else if (i<Banz && B[i] != NULL){
               printf("there does not exists an inhomogenous solution\n");
            }
            if(X[1] != NULL)
            {
               put_mat(X[1], NULL, "homogenous solutions as columns", 0);
               free_mat(X[1]);
            }
            free(X);
        }

        for (i=0;i<Aanz;i++) free_mat(A[i]);
        for (i=0;i<Banz;i++) free_mat(B[i]);
        free(A);
        if (B!=NULL) free(B);
        if (INFO_LEVEL & 8){
           pointer_statistics(0,0);
        }

   exit(0);
}
