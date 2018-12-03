#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <longtools.h>
#include <tools.h>
#include <zass.h>
#include <datei.h>

int SFLAG;

int main(int argc,char *argv[]){

  matrix_TYP **matinv,
              *relator_input,
             **X,
             **out;

  bravais_TYP *G;

  word *relator;

  int anz_erzeuger,
      anz_relatoren,
      i;

  long dim;

  /* ueberpruefen der uebergebenen flags */
  /* juergens standart-function */
  read_header(argc,argv);


  if ((is_option('h') && optionnumber('h') == 0)||(FILEANZ <2)){
     printf("Usage: %s 'file1' 'file2' ['file3'] ['file4'] [-s] [-c] [-d] [-t]\n",argv[0]);
     printf("\n");
     printf("file1: matrix_TYP containing a presentation for G.\n");
     printf("file2: bravais_TYP containing the unimodular group G.\n");
     printf("\n");
     printf("Determines isomorphism type of the cohomology group H^1(G,Q^n/Z^n), or\n");
     printf("in slight abuse of notation (since for finite G they are isomorphic)\n");
     printf("the group of extensions H^2(G,Q^n/Z^n).\n");
     printf("\n");
     printf("Options:\n");
     printf("-s    : writes the system of equations nessecary to\n");
     printf("      : determine the H^1(G,**) to file3\n");
     printf("-c    : writes the cocycles to file3.co or (or file4.co if\n");
     printf("      : used together with -s) and the denomiantors to\n");
     printf("      : file3.de (or file4.de respectively)\n");
     printf("-d    : writes a different output format, which is meant for\n");
     printf("      : developing other programs, so do not us this option.\n");
     printf("-t    : give transformation matrix.\n");
     printf("\n");
     printf("Cf. Extensions/Vectorsystems, Presentation.\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
  }

  if (is_option('h')) { INFO_LEVEL = optionnumber('h'); }
  if (INFO_LEVEL & 8) { SFLAG = 1; }

  /* einlesen der matrizen */
  X = mget_mat(FILENAMES[0],&i);
  if (i>1){
    fprintf(stderr,"you should specify a single matrix for the presentation\n");
    exit(3);
  }
  relator_input = X[0];
  free(X);

  /* read the group */
  G = get_bravais(FILENAMES[1]);
  anz_erzeuger = G->gen_no;

  /* bereitstellen eines pointers fuer die inversen von mat */
  matinv = (matrix_TYP **) calloc(anz_erzeuger,sizeof(matrix_TYP *));

  /* speicher fuer die worte */
  relator = (word *) calloc(relator_input->rows,sizeof(word));

  /* konvertieren der inputmatrix in relator-format */
  for (i=0;i<relator_input->rows;i++){
    matrix_2_word(relator_input,relator+i,i);
  }

  out = cohomology(&dim,G->gen,matinv,relator,anz_erzeuger,relator_input->rows);

  if (dim>0) printf("DIM = %ld\n", dim);

  if (is_option('t') ||
      is_option('d')){
     printf("#3\n");
  }

  if (is_option('d')){
     put_mat(out[0],NULL,"cozykel (columns, for all generators underneath)",2);
     put_mat(out[1],NULL,"devisors",2);
     put_mat(out[2],NULL,"tranformation matrix",2);
     exit(0);
  }

  for (i=out[1]->cols-1;i>=0;i--){
     if (out[1]->array.SZ[i][i] == 1 ||
         out[1]->array.SZ[i][i] == -1 ||
         out[1]->array.SZ[i][i] == 0){
         kill_col(out[1],i);
         kill_row(out[1],i);
     }
  }

  /* now decide whether we really have a cohomology group */
  if (out[0]->cols > 0){
     if (is_option('t')){
        printf("#3\n");
        put_mat(out[0],NULL,
            "cozykel (columns, for all generators underneath)",2);
        put_mat(out[1],NULL,"devisors",2);
     }
     else{
        printf("Isomorphism type of H^1(G,Q^n/Z^n): ");
        printf("C_%d",out[1]->array.SZ[0][0]);
        for (i=1;i<out[0]->cols;i++){
           printf(" X ");
           printf("C_%d",out[1]->array.SZ[i][i]);
        }
        printf("\n");
     }
  }
  else{
     printf("The group has trivial cohomology with this lattice.\n");
     exit(0);
  }

  if (is_option('t')) put_mat(out[2],NULL,"tranformation matrix",2);

  /* freigeben des Speichers fuer die matrizen */
  for (i=0;i<anz_erzeuger;i++){
    if (matinv[i] != NULL){
      free_mat(matinv[i]);
    }
  }
  free(matinv);
  free_mat(out[0]);
  free_mat(out[1]);
  free_mat(out[2]);
  free(out);
  free_bravais(G);

  for (i=0;i<relator_input->rows;i++){
    wordfree(relator+i);
  }
  free(relator);

  free_mat(relator_input);

  /* making some diagnostics */
  if (INFO_LEVEL & 8){
     pointer_statistics(0,0);
  }

  exit(0);
} /* main */

