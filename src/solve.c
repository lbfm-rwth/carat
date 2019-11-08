#error THIS FILE IS NOT USED

#include"typedef.h"

main (int argc, char *argv[])
{

	matrix_TYP **Mat, *X;
        int i, anz;

        read_header(argc, argv);
        if(FILEANZ != 1)
        {
           printf("\n");
           printf("usage:  Solve 'file',\n");
           printf("where 'file' contains a set of matrices.\n");
           printf("\n");
           printf("Solves the linear equation A * X = 0\n");
           printf("for each matrix A in file.\n");
           printf("\n");
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
            X = solve_mat(Mat[i]);
            put_mat(X, NULL, "", 0);
            free_mat(X);
        }

   exit(0);
}
