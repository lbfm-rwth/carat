#include "typedef.h"
#include "matrix.h"
#include "bravais.h"



main (int argc, char *argv[])

{
  bravais_TYP *S;
  bravais_TYP *B;
  matrix_TYP *X;
  int i,j;

  read_header(argc, argv);
  if(FILEANZ != 1)
  {
    printf("usage:   gittstab 'file',\n");
    printf("where 'file' contains a group given as bravais_TYP\n");
    if (is_option('h')){
       exit(0);
    }
    else{
       exit(31);
    }
  }
/***********
  B = (bravais_TYP *) malloc(sizeof(bravais_TYP));
**************/
  B = get_bravais(FILENAMES[0]);
     for(j=0; j<B->zentr_no; j++)
     {
        X = B->zentr[j];
        S = Z_class(B, X);
        put_bravais(S, NULL, "not almost decomposable bravais-group");
        free_bravais(S);
     }


  exit(0);
}
/*{{{}}}*/
