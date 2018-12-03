#include "typedef.h"
#include "getput.h"
#include "bravais.h"

int main (int argc, char *argv[])
{

	matrix_TYP **Mat;
        matrix_TYP **F, *V;
        bravais_TYP *G;
        int i, d, Fanz, anz;

        extern char **FILENAMES;
        extern int FILEANZ;

	extern matrix_TYP **mget_mat();
	extern bravais_TYP *get_bravais();
        extern matrix_TYP *init_mat();
	extern void put_mat ();
        extern void form_to_vec();

        read_header(argc, argv);
        if(FILEANZ != 2 || (is_option('h') && optionnumber('h') == 0))
        {
           printf("Usage: %s 'file1' 'file2' [-m] [-d]\n",argv[0]);
           printf("\n");
           printf(" file1: matrix_TYP\n");
           printf(" file2: matrix_TYP or bravais_TYP.\n");
           printf("\n");
           printf("For each matrix A in file1 a vector V is calculated\n");
           printf("with the following property:\n");
           printf("  A = 1/V[NO+1] * (V[1] * F_1 + V[2] * F_2 +...+ V[NO] * F_NO),\n");
           printf("Where is F_i are the matrices in file2 if file2 is a\n");
           printf("matrix_TYP, otherwise are the matrices describing the\n");
           printf("form space of the bravais group in file2.\n");
           printf("CAUTION: if not used with the option -d, the denominator\n");
           printf("         is not printed, so you will get a vector with\n");
           printf("         only NO columns.\n");
           printf("\n");
           printf("Options:\n");
           printf(" -m: Use a modular (but exact) method to calculate \n");
           printf("     the result. The result is calculated for a couple\n");
           printf("     of big primes and fitted together with the chinese\n");
           printf("     remainder theorem.\n");
           printf("     This method is much faster and avoids trouble with\n");
           printf("     overflow, BUT IS ABLE TO HANDLE THE CASE V[NO+1] == 1\n");
           printf("     ONLY (AND WILL RUN INTO AN INFINITE LOOP OTHERWISE).\n");
           printf(" -d: give the denominator in the NO+1-th colunm.\n");
           printf("\n");
           printf("\n");
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
        V = init_mat(anz, Fanz+1, "");
        if(is_option('m') == TRUE)
        {
          for(i=0;i<anz;i++)
            form_to_vec_modular(V->array.SZ[i], Mat[i], F, Fanz);
          V->cols--;
        }
        else
        {
          for(i=0;i<anz;i++)
          {
            form_to_vec(V->array.SZ[i], Mat[i], F, Fanz, &d);
            V->array.SZ[i][Fanz] = d;
          }
          if(is_option('d'))
            V->cols = Fanz+1;
          else
            V->cols = Fanz;
        }
        put_mat(V, NULL, "matrices written as linear combination", 0);


   exit(0);
}
