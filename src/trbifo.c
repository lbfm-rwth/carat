#include"typedef.h"
#include"getput.h"
#include"bravais.h"


int INFO_LEVEL;

int main (int argc, char *argv[])
{
   bravais_TYP *G, *Gtr;
   matrix_TYP *A, *P, *S;
   int Pmin;

        read_header(argc, argv);
        if(FILEANZ != 2)
        { 
	  printf("Usage: %s 'file1' 'file2'\n",argv[0]);
	  printf("\n");
	  printf("file1: bravais_TYP of the finite unimodular group G.\n");
	  printf("file2: bravais_TYP G^{tr}\n");
	  printf("\n");
	  printf("Calculates the trace bilinear form on the space of invariant forms\n");
	  printf("for G and G^{tr} respectively.\n");
	  printf("If the space of invariant forms is not given for either of the groups, it\n");
	  printf("will be calculated.\n");
	  printf("\n");
	  printf("Cf: Form_elt, Tr_bravais.\n");
          if (is_option('h')){
             exit(0);
          }
          else{
             exit(31);
          }
        }
        G = get_bravais(FILENAMES[0]);
        Gtr = get_bravais(FILENAMES[1]);

	if (G->form_no == 0){
		G->form = formspace(G->gen,G->gen_no,1,&G->form_no);
	}
	if (Gtr->form_no == 0){
		Gtr->form = formspace(Gtr->gen,Gtr->gen_no,1,&Gtr->form_no);
	}

        S = trace_bifo(G->form, Gtr->form, G->form_no);
        put_mat(S, NULL, "tr_bifo", 2);

   exit(0);
}
