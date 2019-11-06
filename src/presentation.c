#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <bravais.h>
#include <datei.h>
#include <longtools.h>
#include <presentation.h>

#define DEBUG FALSE



int main(int argc,char **argv){

  bravais_TYP *G;

  matrix_TYP **base,
              *M;

  bahn **strong;

  int i,
      siz,
      OPT[6];

  char comment[1000];


  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) || (FILEANZ < 1)){
     printf("Usage: %s 'file' [-D]\n",argv[0]);
     printf("\n");
     printf("file: bravais_TYP containing a finite matrix group G.\n");
     printf("\n");
     printf("Calculates a presentation of the finite group G.\n");
     printf("\n");
     printf("Options:\n");
     printf("-D  : meant for debugging. Do not use.\n");
     printf("\n");
     printf("Cf. Is_finite\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }

  INFO_LEVEL = optionnumber('h');

  G = get_bravais(FILENAMES[0]);

  base = get_base(G);

  strong = strong_generators(base,G,TRUE);

  if (DEBUG){
    check_base(strong,G);
  }

  siz = G->order = size(strong);

  if (is_option('D'))
     OPT[0] = 1;
  else
     OPT[0] = 0;

  M = pres(strong,G,OPT);


  if (is_option('D')){
    printf("Size(G) = %d;\n",siz);
    printf("quit;\n");
  }
  else{
    sprintf(comment,"presentation for group in %s",FILENAMES[0]);
    put_mat(M,0,comment,0);
  }

  free_mat(M);
  for (i=0;i<G->dim;i++){
     free_mat(base[i]);
     free_bahn(strong[i]);
     free(strong[i]);
  }
  free(strong);
  free(base);
  free_bravais(G);

  exit(0);

} /* main */

