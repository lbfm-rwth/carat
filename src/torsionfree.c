#include "typedef.h"
#include "matrix.h"
#include "getput.h"
#include "contrib.h"


int INFO_LEVEL;
extern int SFLAG;

void main(int argc,
          char *argv[]){

  bravais_TYP *R;

  int *K,
       order,
       number;

  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h') == 0) || FILEANZ == 0){

     printf("Usage: %s 'file'\n",argv[0]);
     printf("\n");
     printf("file: bravais_TYP containg generators for a space group\n");
     printf("\n");
     printf("Decides whether the given space group is a torsion\n");
     printf("free space group.\n");
     printf("If so, it additionally decides whether it has trivial center.\n");
     printf("\n");
     printf("CAUTION: The program assumes that the translation lattice\n");
     printf("         is in fact Z^n. This is not checked.\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }

  }

  INFO_LEVEL = optionnumber('h');
  if (INFO_LEVEL == 8){
     SFLAG = 1;
  }

  R = get_bravais(FILENAMES[0]);
  K = torsionfree(R,&order,&number);

 
  if (K[0] == TRUE){ 
     printf("The group in %s is torsion free",FILENAMES[0]);
     if (K[1] == TRUE){
        printf(" with trivial center\n");
     }
     else{
        printf("\n");
     }
  }
  else{
     printf("The group in %s is not torsion free\n",FILENAMES[0]);
  }
  printf("The order of the point group is %d,",order);
  printf(" and it has %d conjugacy classes\n",number);

  free_bravais (R);
  free(K);


  if (INFO_LEVEL == 8) pointer_statistics(0,0);

  exit(0);
}
