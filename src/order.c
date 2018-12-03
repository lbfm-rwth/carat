#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <datei.h>

int INFO_LEVEL;
extern int SFLAG;


int main(int argc,char **argv){

  bravais_TYP *G;

  matrix_TYP **base;

  bahn **strong;

  int i,
      j,
      l,
      siz,
      anz;

  char comment[1000];


  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) || (FILEANZ < 1)){
     printf("Usage: %s 'file' [-o]\n",argv[0]);
     printf("\n");
     printf("file: bravais_TYP containing a finite matrix group G.\n");
     printf("\n");
     printf("Calculates the order of G via an Algorithm by Schreier and Sims.\n");
     printf("Options:\n");
     printf("\n");
     printf("-O   : writes the group to stdout again, including the order and\n");
     printf("       a factorisation of it.\n");
     printf("-o   : write a factorisation and the order to stdout. This can\n");
     printf("       be used to append it to a bravais_TYP.\n");
     printf("\n");
     printf("WARNING: THE PROGRAM WILL TERMINATE IFF THE GROUP IS FINITE.\n");
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

  if (INFO_LEVEL & 12){
     SFLAG = 1;
  }

  G = get_bravais(FILENAMES[0]);

  base = get_base(G);

  strong = strong_generators(base,G,FALSE);

  siz = G->order = size(strong);

  memset(G->divisors,0,100*sizeof(int));

  for (i=2;i<100;i++){
     if ((siz % i) == 0){
       siz /= i;
       G->divisors[i]++;
       i--;
     }
  }

  if (is_option('o')){
     fput_order(stdout,G->divisors,G->order);
  }
  if (is_option('O')){
     put_bravais(G,NULL,NULL);
  }
  fprintf(stderr,"The order of the group is %d\n",G->order);

  free_bravais(G);
  for (i=0;i<G->dim;i++){
     free_mat(base[i]);
     free_bahn(strong[i]);
     free(strong[i]);
  }
  free(strong);
  free(base);

  if (INFO_LEVEL & 12){
     pointer_statistics(0,0);
  }

  exit(0);

} /* main */

