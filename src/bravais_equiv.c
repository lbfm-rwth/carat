#include "typedef.h"
#include "longtools.h"
#include "getput.h"
#include "bravais.h"
#include "symm.h"
#include "autgrp.h"
#include "matrix.h"
#include "voronoi.h"
#include "polyeder.h"
#include "datei.h"

int INFO_LEVEL;
extern int SFLAG;

int main (int argc, char *argv[])
{

  bravais_TYP *G,
              *G_tr,
              *H,
              *H_tr;

  matrix_TYP *erg;

  char comment[1000];

  int i,
      tmp;

  read_header(argc, argv);
  if ((FILEANZ < 2) || (is_option('h') && optionnumber('h') ==0)){
    printf("Usage: %s 'file1' 'file2'\n",argv[0]);
    printf("\n");
    printf("file1: bravais_TYP containing G.\n");
    printf("file2: bravais_TYP containing H.\n");
    printf("\n");
    printf("Tests whether the BRAVAIS GROUPS of the groups\n");
    printf("G and H respectively are conjugated in GL_n(Z).\n");
    printf("If so, it returns a conjugating matrix X which conjugates\n");
    printf("the BRAVAIS GROUPS, ie. X^1 B(G) X = B(H).\n");
    printf("\n");
    printf("WARNING: The procedure may involve calculating the normalizer\n");
    printf("         of the groups. This may be very time comsuming,\n");
    printf("         especially when both groups are <-I_n>, where n>5.\n");
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

  /* we might not deal with th whole bravais group */
  G->order = H->order = 0;

  /* paranoia setting: recalculate the formspace because it has to be a
     Z-basis */
  if (G->form !=NULL){
     for (i=0;i<G->form_no;i++){
        free_mat(G->form[i]);
     }
     free(G->form);
  }
  G->form = formspace(G->gen,G->gen_no,1,&tmp);
  G->form_no = tmp;
  if (H->form !=NULL){
     for (i=0;i<H->form_no;i++){
        free_mat(H->form[i]);
     }
     free(H->form);
  }
  H->form = formspace(H->gen,H->gen_no,1,&tmp);
  H->form_no = tmp;

  G_tr = tr_bravais(G,1,FALSE);
  H_tr = tr_bravais(H,1,FALSE);

  /* output for debugging purposes
  put_bravais(G,NULL,NULL);
  put_bravais(H,NULL,NULL); */

  erg = is_z_equivalent(G,G_tr,H,H_tr);

  if (erg == NULL){
     printf("the bravais groups are not conjugated\n");
  }
  else{
     sprintf(comment,"conjugates the group of %s in the group of %s",
                      FILENAMES[0],FILENAMES[1]);
     put_mat(erg,NULL,comment,2);
     free_mat(erg);
  }

  /* output the groups again, just to make sure we didn't change
     them */
  if (INFO_LEVEL == 5){
     put_bravais(H,NULL,NULL);
     put_bravais(H_tr,NULL,NULL);
  }

  free_bravais(G);
  free_bravais(G_tr);
  free_bravais(H);
  free_bravais(H_tr);


  if (INFO_LEVEL & 12){
     pointer_statistics(0,0);
  }
  exit(0);
}
