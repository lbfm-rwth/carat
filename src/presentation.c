#include "typedef.h"
#include "matrix.h"
#include "symm.h"
#include "getput.h"
#include "bravais.h"
#include "sort.h"
#include "polyeder.h"
#include "presentation.h"
#include "tools.h"
#include "tietzetrans.h"
#include "datei.h"

extern int SFLAG;

void main (int argc, char *argv[])
{
  bravais_TYP 		*G;
  presentation_TYP	*erg;

  read_header(argc, argv);
  if(FILEANZ != 1)
  {
    printf("\n");
    printf("Usage\n");
    printf("  Presentation file [-G]\n");
    printf(" \n");
    printf(" where file contains a bravais_TYP describing\n");
    printf(" a FINITE INTEGRAL group.\n");
    printf(" \n");
    printf(" Calculates a presentation of the group in file, and writes\n");
    printf(" it to the stdout.\n");
    printf(" \n");
    printf(" The options are:\n");
    printf("\n");
    printf(" -G   : Output the presentation in Gap format on a file\n");
    printf("        called <file>.gap.short .\n");
    printf(" -h   : Gives you this help.\n");
    printf("\n");
    printf("\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }

  /* setting SFLAG according to optionnumber('h') */
  if (is_option('h') && optionnumber('h') == 8){
     SFLAG = 1;
  }

  G = get_bravais(FILENAMES[0]);

  erg = presentation_point_grp(G);

  put_presentation(erg,NULL,"C");

  free_presentation(erg);

  free_bravais(G);

  if (SFLAG == 1){
     pointer_statistics(0,0);
  }

  exit(0);
}
