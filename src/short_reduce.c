#include"typedef.h"
#include"getput.h"
#include"matrix.h"
#include"symm.h"
#include"reduction.h"




int main (int argc, char *argv[])
{

	matrix_TYP **Mat, *T, *Red, *SV;
        int i, anz;
        int min;

        extern char **FILENAMES;
        extern int FILEANZ;

	extern matrix_TYP *short_reduce();

        read_header(argc, argv);
        if(FILEANZ != 1)
        {
           printf("usage:  short_reduce 'file',\n");
           printf("where 'file' contains a set of positiv semidefinite matrices.\n");
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
            T = init_mat(Mat[i]->cols, Mat[i]->cols, "1");
            SV = shortest(Mat[i], &min);
            Red = short_reduce(Mat[i], SV, T);
            put_mat(Red, NULL, "reduced matrix", 2);
            if(is_option('t') == TRUE)
             put_mat(T, NULL, "Transformation matrix", 0);
            free_mat(Red);
            free_mat(T);
            free_mat(SV);
        }

   exit(0);
}
