#include "typedef.h"
#include "matrix.h"
#include "getput.h"
#include "contrib.h"
#include "datei.h"

int INFO_LEVEL;
extern int SFLAG;

int main(int argc,
          char *argv[]){

  bravais_TYP *R;

  int *K,
       order,
       number;

  read_header(argc,argv);

  if ((is_option('h') && optionnumber('h') == 0) || FILEANZ == 0){

     printf("Usage: %s 'file' [-t]\n",argv[0]);
     printf("\n");
     printf("file: bravais_TYP containg generators for a space group\n");
     printf("\n");
     printf("Decides whether the given space group is a torsion\n");
     printf("free space group.\n");
     printf("If so, it additionally decides whether it has trivial center.\n");
     printf("\n");
     printf("CAUTION: The program assumes that the translation lattice\n");
     printf("         is in fact Z^n. This is not checked.\n");
     printf("\n");
     printf("Options:\n");  // Oliver: 11.04.2002
     printf("-t    : Output for TeX. Do not use!\n");
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

  if (!is_option('t')){
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
  }
  else{		// Oliver: 11.04.2002
     if (K[0] == TRUE){
        if (R->dim != 5)
           printf("$\\times$");
        if (K[1] == TRUE)
           printf("$\\times$ ");
     }
/*     else{
        printf(" & ");
     } */
     printf("\n");
  }

  free_bravais (R);
  free(K);


  if (INFO_LEVEL == 8) pointer_statistics(0,0);

  exit(0);
}
