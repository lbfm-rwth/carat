#error THIS FILE IS NOT USED

#include"typedef.h"


main (int argc, char *argv[])
{

	matrix_TYP **x1, **x2, **S, *X;
        int anz;

        read_header(argc, argv);
        if(FILEANZ != 3)
        {
           printf("usage:  hypstab 'file1' 'file2' 'file3' ,\n");
           printf("where 'file1' contains an 1 by n matrix x.\n");
           printf("where 'file2' contains an 1 by n matrix y.\n");
           printf("and   'file3' contains a n by n matrix S.\n");
           printf("of sylvester Typ (-1, n-1) s.t. xSx^{tr} < 0 and ySy^{tr} < 0\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }
        x1 = mget_mat(FILENAMES[0], &anz);
        x2 = mget_mat(FILENAMES[1], &anz);
        S = mget_mat(FILENAMES[2], &anz);
        X = hyperbolic_isometry(x1[0], x2[0], S[0]);
        if(X == NULL)
         printf("There exists no isometry\n");
        else
           put_mat(X, NULL, "isometry", 0);


   exit(3);
}
