#include "typedef.h"
#include "getput.h"
#include "matrix.h"


int main (int argc, char *argv[])
{

	matrix_TYP *Mat, *R;
        matrix_TYP *SV;
        int denomin, length, anz, i;

        extern char **FILENAMES;
        extern int FILEANZ;

	extern matrix_TYP *get_mat ();
	extern matrix_TYP *rest_short();
	extern void put_mat ();

        read_header(argc, argv);
        if(FILEANZ != 2)
        {
          printf("usage: Rest_short 'file1' 'file2' [-l=n -m=n -f -c -n],\n");
          printf("where 'file1' contains a positive definite symmetric matrix\n");
          printf("and   'file2' contains a matrix\n");
          if (is_option('h')){
             exit(0);
          }
          else{
             exit(31);
          }
        }
	Mat = get_mat (FILENAMES[0]);
        R = get_mat(FILENAMES[1]);
        length = optionnumber('l');
        if(length == 0)
          length = 1;
        denomin = optionnumber('d');
        if(denomin == 0)
           denomin = 1;
        for(i=0;i<R->rows;i++)
        {
          SV = rest_short(Mat, R->array.SZ[i], R->kgv, length, denomin,is_option('f'), is_option('c'), &anz);
          if(is_option('c') == TRUE)
          {
            printf("#%d  Anzahl der Vektoren mit Quadratabstand kleiner %d/%d\n", anz, denomin, length);
          }
          else
          {
            if(is_option('n') == TRUE)
               SV->cols++;
            put_mat(SV, NULL, "shortest vectors", 0);
            free_mat(SV);
          }
        }

   exit(0);
}
