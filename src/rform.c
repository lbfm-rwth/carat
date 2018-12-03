#include"typedef.h"
#include"matrix.h"
#include"bravais.h"
#include"getput.h"


int main (int argc, char *argv[])
{

	matrix_TYP **Mat, *V;
        bravais_TYP *G;
        int epsilon, Mat_anz;
        int i,j;


        read_header(argc, argv);
        if(FILEANZ != 2 && FILEANZ != 1)
	{
	  printf("Usage:  %s 'file1' ['file2'] [-e=epsilon]\n",argv[0]);
	  printf("\n");
	  printf("file1: bravais_TYP of the finite unimodular group G.\n");
	  printf("file2: (OPTIONAL) matrix_TYP of a single matrix A. If file2 is not\n");
	  printf("       given, A=I_n the identity matrix is assumed.\n");
	  printf("\n");
	  printf("The program calculates the sum of g^{tr}Ag with g in G in a seminumerical\n");
	  printf("way, and divides the results by the greatest common divisor of it's\n");
	  printf("entries.\n");
	  printf("Note that this gives a way to construct a positive definite, G-invariant\n");
	  printf("matrix.\n");
	  printf("\n");
	  printf("Options:\n");
	  printf("-e=epsilon : epsilon is a control parameter, controlling the accuracy\n");
	  printf("             of the convergence. It is given as an integer epsilon>=100, and\n");
	  printf("             bigger value stand for a more accurate  convergence.\n");
	  printf("\n");
	  printf("Cf. Invar_space\n");
          if (is_option('h')){
             exit(0);
          }
          else{
             exit(31);
          }
	}
        G = get_bravais(FILENAMES[0]);
        if(FILEANZ == 2)
	  Mat = mget_mat (FILENAMES[1], &Mat_anz);
        else
        {
          if( (Mat = (matrix_TYP **)malloc(1 *sizeof(matrix_TYP *))) == NULL){
            printf("malloc of 'Mat' in main program failed\n");
            exit(2);
          }
          Mat[0] = init_mat(G->dim, G->dim, "1");
          Mat_anz  = 1;
        }
        epsilon = optionnumber('e');
        if(epsilon == 0)
         epsilon = 100;
        if(Mat_anz != 1)
         printf("#%d\n", Mat_anz);
        for(i=0;i<Mat_anz;i++)
        {
           V = rform(G->gen, G->gen_no, Mat[i], epsilon);
             put_mat(V, NULL, "", 2);
           free_mat(V);
        }

   exit(0);
}
