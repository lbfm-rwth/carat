#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <datei.h>

int INFO_LEVEL;
extern int SFLAG;

int main(int argc,char **argv){

  int i,
      tmp;

  bravais_TYP *G;

  matrix_TYP **base = NULL;

  bahn **strong = NULL;

  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h')==0) || (FILEANZ < 1)){
     printf("Usage: %s 'file'\n",argv[0]);
     printf("\n");
     printf("file: bravais_TYP containing the FINITE group G.\n");
     printf("\n");
     printf("Tries to reduce the number of generators for the group in file. It\n");
     printf("will output a bravais_TYP, where the generators are a subset of the\n");
     printf("given generator such that they generate the same group.\n");
     printf("\n");
     printf("Cf. Is_finite, Order.\n");
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
  memset(G->divisors,0,100);

  base = get_base(G);

  tmp = G->order = red_gen(G,base, &strong, 0);

  for (i=2;i<101;i++){
     if (tmp % i == 0){
        tmp = tmp/i;
        G->divisors[i]++;
        i--;
     }
  }
  put_bravais(G,NULL,NULL);

  for (i=0;i<G->dim;i++){
     free_bahn(strong[i]);
     free_mat(base[i]);
     free(strong[i]);
  }
  free(strong);
  free(base);
  free_bravais(G);

  if (INFO_LEVEL & 12){
     pointer_statistics(0,0);
  }

  exit(0);
} /* main */
