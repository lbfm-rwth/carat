#include "typedef.h"
#include "getput.h"
#include "matrix.h"

int main (int argc, char *argv[])
{

	matrix_TYP **Mat, **V;
        bravais_TYP *G;
        int V_anz, Mat_anz;
        int i,j;

        extern char **FILENAMES;
        extern int FILEANZ;

	extern matrix_TYP **mget_mat();
	extern bravais_TYP *get_bravais();
        extern matrix_TYP **voronoi_vertices();
	extern void put_mat ();

        read_header(argc, argv);
        if(FILEANZ != 2)
	{
          printf("usage:  Vor_vertices 'file1' 'file2',\n");
          printf("where 'file1' containes a G-invariant positive definite, symmetric matrix and\n");
          printf("      'file2' containes a group G, given as bravais_TYP\n");
          if (is_option('h')){
             exit(0);
          }
          else{
             exit(31);
          }
	}
	Mat = mget_mat (FILENAMES[0], &Mat_anz);
        G = get_bravais(FILENAMES[1]);
        for(i=0;i<Mat_anz;i++)
        {
           V = voronoi_vertices(Mat[i], G, &V_anz, &j);
           printf("#%d\n", V_anz);
           for(j=0;j<V_anz;j++)
           {
             put_mat(V[j], NULL, "", 0);
             free_mat(V[j]);
           }
           free(V);
        }

   exit(0);
}
