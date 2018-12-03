#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <datei.h>

int INFO_LEVEL;
extern int SFLAG;

matrix_TYP *conjugated(bravais_TYP *G,bravais_TYP *H,
                       matrix_TYP **N,int Nanz,bahn **strong);


int main(int argc,char **argv){

  bravais_TYP *G,
              *H,
              *N;

  matrix_TYP **base,
             *erg;

  bahn **strong;

  int i,
      j,
      l,
      anz;

  char comment[1000];


  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) || (FILEANZ < 3)){
     printf("Usage: Conjugated file1 file2 file3\n");
     printf("\n");
     printf("file1: bravais_TYP containing the finite group G.\n");
     printf("file2: bravais_TYP containing the finite group H.\n");
     printf("file2: bravais_TYP containing the group N.\n");
     printf("\n");
     printf("The program checks whether there is an x in N with the\n");
     printf("property that x H x^{-1} <= G, and returns this x if it\n");
     printf("exists.\n");
     printf("\n");
     printf("WARNING: The program assumes that the orbit of H under N\n");
     printf("         is finite.\n");
     printf("\n");
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
  H = get_bravais(FILENAMES[1]);
  N = get_bravais(FILENAMES[2]);

  base = get_base(G);
  strong = strong_generators(base,G,FALSE);

  erg = conjugated(G,H,N->gen,N->gen_no,strong);

  /* output */
  if (erg != NULL){
     sprintf(comment,"conjugated the group of %s INTO the group of %s",
                      FILENAMES[1],FILENAMES[0]);
     put_mat(erg,NULL,comment,2);
     free_mat(erg);
  }
  else{
     printf("the groups are not conjugated under the group of %s\n",
             FILENAMES[2]);
  }


  free_bravais(G);
  free_bravais(H);
  free_bravais(N);
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

