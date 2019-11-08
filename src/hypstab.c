#error THIS FILE IS NOT USED

#include"typedef.h"


main (int argc, char *argv[])
{

	matrix_TYP **x, **S;
        bravais_TYP *G;
        int anz;

        read_header(argc, argv);
        if(FILEANZ != 2)
        {
           printf("usage:  hypstab 'file1' 'file2 ,\n");
           printf("where 'file1' contains an 1 by n matrix X.\n");
           printf("and   'file2' contains a set of n by n matrix S.\n");
           printf("of sylvester Typ (-1, n-1) s.t. xSx^{tr} < 0\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }
        x = mget_mat(FILENAMES[0], &anz);
        S = mget_mat(FILENAMES[1], &anz);
        G = hyperbolic_stabilizer(x[0], S[0]);
        put_bravais(G, NULL, "hyperbolic stabilizer");


   exit(0);
}
