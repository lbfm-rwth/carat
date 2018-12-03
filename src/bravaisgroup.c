#include "typedef.h"
#include "getput.h"
#include "bravais.h"
#include "symm.h"
#include "datei.h"

int INFO_LEVEL;
extern int SFLAG;

int main (int argc, char *argv[])
{

  bravais_TYP *G,
              *H;

  char comment[1000];

  read_header(argc, argv);
  if ((FILEANZ != 1) || (is_option('h') && optionnumber('h') ==0)){
    printf("Usage: %s 'file'\n",argv[0]);
    printf("\n");
    printf("file: bravais_TYP containing the finite unimodular group G.\n");
    printf("\n");
    printf("Calculates generators of the Bravais group B(G) of G. If \n");
    printf("the space of G-invariant quadratic forms is given in 'file',\n");
    printf("it relies on its correctness and echoes it in the \n");
    printf("output. Otherwise it is calculated as well.\n");
    printf("\n");
    printf("Cf. Aut_grp (for forms without a group in file).\n");
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

  H = bravais_group(G,FALSE);

  sprintf(comment,"bravais group to %s",FILENAMES[0]);
  put_bravais(H,NULL,comment);

  free_bravais(G);
  free_bravais(H);

  if (INFO_LEVEL & 12){
     pointer_statistics(0,0);
  }
  exit(0);
}
