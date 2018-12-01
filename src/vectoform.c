#include"typedef.h"
#include"getput.h"
#include"matrix.h"
#include"bravais.h"


int main (int argc, char *argv[])
{

	matrix_TYP **Mat;
        matrix_TYP **F, *V;
        bravais_TYP *G;
        int i, d, Fanz, anz;
        char comment[1024];

        read_header(argc, argv);
        if(FILEANZ != 2 || (is_option('h') && optionnumber('h') == 0))
        {
           printf("Usage: %s 'file1' 'file2' \n",argv[0]);
           printf("\n");
           printf("file1: matrix_TYP of 1 by n matrices.\n");
           printf("file2: matrix_TYP or bravais_TYP. If file2 is a matrix_TYP, it has to\n");
           printf("       contain exactly n matrices, if it contains a bravais_TYP its\n");
           printf("       space of invariant forms should be of dimension n.\n");
           printf("\n");
           printf("For each row vector V in file1 a matrix A is calculated\n");
           printf("with the following property:\n");
           printf("  A =  (V[1] F_1 + V[2] * F_2 +...+ V[NO] * F_NO),\n");
           printf("Where is F_i are the matrices in file2 if file2 is a\n");
           printf("matrix_TYP, otherwise are the matrices describing the\n");
           printf("space of invariant forms of the bravais group in file2.\n");
           printf("\n");
           printf("Cf. Formtovec\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }
	Mat = mget_mat (FILENAMES[0], &anz);
        G = get_bravais(FILENAMES[1]);
        if(G->form_no > 0)
           { F = G->form; Fanz = G->form_no;}
        else
           { F = G->gen;  Fanz = G->gen_no;}

        for (i=0;i<anz;i++){
           if (Mat[i]->cols < Fanz){
              fprintf(stderr,"uncompatible dimensions\n");
              exit(3);
           }
           if (!Mat[i]->flags.Integral){
              fprintf(stderr,"expecting an integral matrix\n");
              exit(3);
           }
           if (Mat[i]->rows != 1){
              fprintf(stderr,"expecting a matrix with a single row\n");
              exit(3);
           }
           V = vec_to_form(Mat[i]->array.SZ[0],F,Fanz);
           sprintf(comment,"linear combination of matrices in %s",FILENAMES[1]);
           put_mat(V,NULL,comment,2);
           free_mat(V);
        }

  exit(0);
}
