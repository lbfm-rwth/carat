#include "typedef.h"
#include "getput.h"
#include "matrix.h"
#include "symm.h"


int main (int argc, char *argv[])
{

	matrix_TYP **Mat, *S;
        int i,j,n, anz;
        int test;

        extern char **FILENAMES;
	extern int FILEANZ;

	extern matrix_TYP **mget_mat ();
        extern matrix_TYP *dsylv();
	extern void put_mat ();

        read_header(argc, argv);
        if(FILEANZ != 1)
        {
           printf("Usage: %s 'file' [-t]\n",argv[0]);
           printf("\n");
           printf("file:  matrix_TYP containing a set of symmetric matrices\n");
           printf("\n");
           printf("Calculates the signature (Sylvester type) of the matrices in file.\n");
           printf("For a matrix congruent to diag(1,...,1,-1,...,-1,0,...0) over the\n");
           printf("reals it returns the appropriate number of 1's, -1's, and 0's.\n");
           printf("\n");
           printf("Options:\n");
           printf("-t:   it decides only which matrices are positive or negative (semi-)definite.\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }
	Mat = mget_mat (FILENAMES[0], &anz);
        if(anz > 1)
          printf("#%d\n", anz);
        for(i=0;i<anz;i++)
        {
            n = Mat[i]->cols;
            Check_mat(Mat[i]);
            if(Mat[i]->flags.Symmetric == FALSE)
              printf("matrix %d is not symmetric\n", i+1);
            else
            {
              if(is_option('t') == TRUE)
              {
                 test = definite_test(Mat[i]);
                 if(test == 2)
                 printf("matrix no. %d is positiv definite\n", (i+1));
                 if(test == 1)
                 printf("matrix no. %d is positiv semidefinite, but not definite\n", (i+1));
                 if(test == -2)
                 printf("matrix no. %d is negativ definite\n", (i+1));
                 if(test == -1)
                 printf("matrix no. %d is negativ semidefinite, but not definite\n", (i+1));
                 if(test == -3)
                 printf("matrix no. %d is indefinite\n", (i+1));
                 if(test == 0)
                 printf("matrix no. %d is zero\n", (i+1));
              }
              else
              {
                S = dsylv(Mat[i]);
                put_mat(S, NULL, "", 2);
                free_mat(S);
              }
            }
        }

   exit(0);
}
