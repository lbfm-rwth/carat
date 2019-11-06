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


int main (int argc, char *argv[])
{
  bravais_TYP 		*G;
  matrix_TYP		*vec;
  matrix_TYP		*Form;
  polyeder_TYP		*Pol;

        read_header(argc, argv);
        if(FILEANZ != 3)
        {
   printf("\n");
   printf("Usage\n");
   printf("  %s file1 file2 file3 \n", argv[0]);
   printf(" \n");
   printf(" where file1 contains a bravais_TYP describing\n");
   printf(" a set of affine matrices generating a space group.\n");
   printf(" \n");
   printf(" where file2 contains a matrix_TYP describing\n");
   printf(" an affine vector which is a starting point for the algorithm\n");
   printf(" \n");
   printf(" where file3 contains a matrix_TYP describing a positive definite,\n");
   printf(" invariant form for R \n");
   printf(" \n");
   printf(" Calculates a fundamental polyhedron for the group in file1,\n");
   printf(" and writes it to the stdout.\n");
 
   printf(" \n");
   printf(" The options are:\n");
   printf("\n");
   printf(" -h   : Gives you this help.\n");
   printf("\n");

          if (is_option('h')){
             exit(0);
          }
          else{
             exit(31);
          }
        }

        G = get_bravais(FILENAMES[0]);
        vec = get_mat(FILENAMES[1]);
        Form = get_mat(FILENAMES[2]);

        Pol = fub(vec, G, Form);

        put_polyeder(Pol);
        free_polyeder(Pol);

        free_bravais(G);
        free_mat(vec); vec = NULL;
        free_mat(Form); Form = NULL;

   printf("\n");
}
