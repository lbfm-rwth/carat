#include "typedef.h"
#include "getput.h"
#include "polyeder.h"

int INFO_LEVEL;



int main (int argc, char *argv[])
{

	matrix_TYP *Mat;
        polyeder_TYP *P;
        wall_TYP **w;
        int i,j, anz;

        extern char **FILENAMES;
        extern int FILEANZ;

	extern matrix_TYP *get_mat ();
        extern polyeder_TYP *first_polyeder();
        extern wall_TYP *init_wall();

        read_header(argc, argv);
        if(FILEANZ != 1)
        {
          printf("usage: Polyeder 'file',\n");
          printf("where 'file' contains a matrix, whose rows are interpreted as linear inequalities\n");
          if (is_option('h')){
             exit(0);
          }
          else{
             exit(31);
          }
        }
	Mat = get_mat (FILENAMES[0]);
        w = (wall_TYP **)malloc(Mat->rows *sizeof(wall_TYP *));
        for(i=0;i<Mat->rows;i++)
        {
            w[i] = init_wall(Mat->cols);
            for(j=0;j<Mat->cols;j++)
              w[i]->gl[j] = Mat->array.SZ[i][j];
        }
        P = first_polyeder(w, Mat->rows);
        for(i=0;i<Mat->rows;i++)
        {
           j =  refine_polyeder(P, w[i]);
        }
        put_polyeder(P);

   exit(0);
}
