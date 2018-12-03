#include"typedef.h"
#include"symm.h"
#include"getput.h"


int main (int argc, char *argv[])
{

	matrix_TYP *Mat;
        matrix_TYP *SV;
        int min;

        read_header(argc, argv);
        if(FILEANZ != 1)
        {
           printf("Usage: %s 'file' [-n]\n",argv[0]);
           printf("\n");
           printf("file: matrix_TYP containing a positive definite symmetric matrix.\n");
           printf("\n");
           printf("Calculates the shortest vectors of the given form. They are output as\n");
           printf("a matrix whose rows are the shortest vectors of the form.\n");
           printf("\n");
           printf("Options:\n");
           printf("-n     : Echo the norm of the vectors in the last column of the vector.\n");
           printf("\n");
           printf("Cf: Isometry, Aut_grp, Short.\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }
	Mat = get_mat (FILENAMES[0]);
        SV = shortest(Mat, &min);
        if(is_option('n') == TRUE)
           SV->cols++;
        put_mat(SV, NULL, "shortest vectors", 0);

   exit(0);
}
