#include"typedef.h"
#include"tools.h"
#include"matrix.h"
#include"getput.h"
#include"longtools.h"




int main (int argc, char *argv[])
{
  matrix_TYP **M, *Trf, *E;
  int Manz, i;

  extern char **FILENAMES;
  extern int FILEANZ;


  read_header(argc, argv);
  if(FILEANZ != 1)
  {
     printf("Usage: %s 'file' [-t]\n",argv[0]);
     printf("\n");
     printf("file:  matrix_TYP  containing a set of matrices\n");
     printf("\n");
     printf("Calculates the Smith normal form of the matrices\n");
     printf("given in file, resp. the isomorphism type of the factor group of\n");
     printf("two lattices, where the matrix is interpreted as the coordinate\n");
     printf("columns of a generating set of the second lattice in terms of a\n");
     printf("basis of the first one.\n");
     printf("\n");
     printf("CAUTION: If the given matrix is rational, the\n");
     printf("         Smith normal form is calculated for the least\n");
     printf("         integral multiple of it.\n");
     printf("\n");
     printf("Options:\n");
     printf("-t   : gives the left tranformation matrix as well.\n");
     printf("\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }

  M = mget_mat(FILENAMES[0], &Manz);
  for(i=0;i<Manz;i++)
  {
    if (!M[i]->flags.Integral) rat2kgv(M[i]);

    if(is_option('t'))
      Trf = init_mat(M[i]->rows, M[i]->rows, "1");
    else
      Trf = NULL;
    E = long_elt_mat(Trf, M[i], NULL);
    put_mat(E, NULL, "ELT of matrix", 0);
    free_mat(E);
    if(is_option('t'))
    {
      put_mat(Trf, NULL, "Transformtion matrix", 0);
      free_mat(Trf);
    }
  }


  exit(0);
}
