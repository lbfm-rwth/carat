#include"typedef.h"
#include"getput.h"
#include"symm.h"


int main (int argc, char *argv[])
{

	matrix_TYP *Mat;
        matrix_TYP *SV;
        int min, length, anz;

        read_header(argc, argv);
        if(FILEANZ != 1)
        {
           printf("Usage: %s 'file' [-l=n] [-m=n] [-f] [-c] [-n]\n",argv[0]);
           printf("\n");
           printf("file: matrix_TYP containing a positive definite symmetric matrix.\n");
           printf("\n");
           printf("Calculates the short vectors of the given form. They are output as\n");
           printf("a matrix whose rows are the short vectors of the form.\n");
           printf("\n");
           printf("Options:\n");
           printf("-l=n   : Searches for vectors up to the given norm n.\n");
           printf("-m=n   : Searches for vectors with minimal norm n.\n");
           printf("-f     : Stop immediately if a vector with the given property is found.\n");
           printf("-c     : Just echo the number of vectors with the given property.\n");
           printf("-n     : Echo the norm of the vectors in the last column of the vector.\n");
           printf("\n");
           printf("Cf: Isometry, Aut_grp, Shortest.\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }
	Mat = get_mat (FILENAMES[0]);
        length = optionnumber('l');
        if(length == 0)
          length = 1;
        min = optionnumber('m');
        SV = short_vectors(Mat, length, min,is_option('f'), is_option('c'), &anz);
        if(is_option('c') == TRUE)
        {
          printf("#%d  Anzahl der Vektoren mit Norm zwischen %d und %d\n", anz, min, length);
        }
        else
        {
          if(is_option('n') == TRUE)
             SV->cols++;
          put_mat(SV, NULL, "shortest vectors", 0);
        }

   exit(0);
}
