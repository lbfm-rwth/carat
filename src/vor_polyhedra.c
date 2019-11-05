#include"typedef.h"


main (int argc, char *argv[])
{

	matrix_TYP **Mat;
        bravais_TYP *G, *Gtr;
        polyeder_TYP *pol;
        matrix_TYP *bifo;
        int Mat_anz;

        read_header(argc, argv);
        if(FILEANZ != 3)
	{
          printf("usage:  vor_polyhedra 'file1' 'file2' 'file3',\n");
          printf("where 'file1' containes a group G, given as bravais_TYP\n");
          printf("      'file2' containes the transposed group G^{tr}, given as bravais_TYP\n");
          printf("      'file3' containes a G-perfect matrix \n");
          if (is_option('h')){
             exit(0);
          }
          else{
             exit(31);
          }
	}
        G = get_bravais(FILENAMES[0]);
        Gtr = get_bravais(FILENAMES[1]);
	Mat = mget_mat (FILENAMES[2], &Mat_anz);
        bifo = trace_bifo(G->form, Gtr->form, G->form_no);
        pol = vor_polyhedra(Mat[0], G, Gtr->form, bifo);
        put_polyeder(pol);


   exit(0);
}
