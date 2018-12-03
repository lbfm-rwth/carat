#include "typedef.h"
#include "longtools.h"
#include "getput.h"
#include "bravais.h"
#include "matrix.h"
#include "datei.h"

int INFO_LEVEL;

int main (int argc, char *argv[])
{

  bravais_TYP *G,
              *G_tr;

  matrix_TYP *tmp,
             *tmp2;

  int i,
      j,
      prime=1949,
      eps=100;

  char comment[1000];

  read_header(argc, argv);
  if(FILEANZ < 1)
  {
    printf("Usage: %s 'file' [-i] [-f]\n",argv[0]);
    printf("\n");
    printf("file: bravais_TYP containing the finite unimodular group G.\n");
    printf("\n");
    printf("Calculates the elementatry divisors of the trace bifo\n");
    printf("of the group G given in file, i.e. takes a Z-basis B (resp. B')\n");
    printf("of the invariant forms of G (respectively G^tr), and\n");
    printf("gives the Smith normal form of Trace(B_i * B'_j)_i,j.\n");
    printf("\n");
    printf("Options:\n");
    printf("\n");
    printf("-i  : Use the function invar_space to calculate the space of invarinat forms.\n");
    printf("-f  : Force the formspace of G to be computed again, even if\n");
    printf("      it is given.\n");
    printf("\n");

    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  G = get_bravais(FILENAMES[0]);
  G_tr = (bravais_TYP *) malloc(1 * sizeof(bravais_TYP));
  G_tr->gen = (matrix_TYP **) malloc(G->gen_no*sizeof(matrix_TYP*));

  for (i=0;i<G->gen_no;i++){
     G_tr->gen[i] = tr_pose(G->gen[i]);
  }

  /* calculate the formspace if it's not known or the secure option -f is
     given */
  if ((G->form==NULL) || is_option('f')){
     /* decide which way the user wants to formspace to be computed */
     if (is_option('i')){
        G->form = p_formspace(G->gen,G->gen_no,prime,1,&G->form_no);
        G->form = invar_space(G->gen,G->gen_no,G->form_no,1,eps,&G->form_no);
     }
     else{
        G->form = formspace(G->gen,G->gen_no,1,&G->form_no);
     }
  }

  if (is_option('i')){
     G_tr->form = invar_space(G_tr->gen,G->gen_no,G->form_no+1,1,eps,
                                                    &G_tr->form_no);
  }
  else{
     G_tr->form = formspace(G_tr->gen,G->gen_no,1,&G_tr->form_no);
  }

  tmp = trace_bifo(G->form,G_tr->form,G->form_no);

  tmp2 = long_elt_mat(NULL,tmp,NULL);

  sprintf(comment,"elementary divisors of the trace_bifo for G in %s",
                                                       FILENAMES[0]);
  put_mat(tmp2,NULL,comment,2);

  free_mat(tmp);
  free_mat(tmp2);
  free_bravais(G);


  exit(0);
}
