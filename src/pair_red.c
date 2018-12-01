#include"typedef.h"
#include"getput.h"
#include"matrix.h"
#include"reduction.h"

int main (int argc, char *argv[])
{

	matrix_TYP **Mat, *T, *Red;
        int i, anz;

        read_header(argc, argv);
        if(FILEANZ != 1)
        {
          printf("Usage: %s 'file' [-t]\n",argv[0]);
          printf("\n");
          printf("file: matrix_TYP containing positive definite symmetric matrices.\n");
          printf("\n");
          printf("The program applies the pair reduction algorithm to each matrix in file.\n");
          printf("\n");
          printf("Options:\n");
          printf("-t    : Give the transforming matrices as well, i.e. matrices\n");
          printf("        T satisfying T^tr * F_old * T = F_new.\n");
          printf("\n");
          printf("Cf. Mink_red, Rform, Conj_bravais.\n");
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
            T = init_mat(Mat[i]->cols, Mat[i]->cols, "");
            Red = pair_red(Mat[i], T);
            put_mat(Red, NULL, "reduced matrix", 2);
            if(is_option('t') == TRUE)
             put_mat(T, NULL, "Transformation matrix", 0);
            free_mat(Red);
            free_mat(T);
        }

   exit(0);
}
