#include"typedef.h"
#include"getput.h"
#include"bravais.h"

#include"longtools.h"
#include"symm.h"
#include"base.h"
#include"autgrp.h"
#include"matrix.h"
#include"voronoi.h"
#include"polyeder.h"
#include "datei.h"


int main (int argc, char *argv[])
{

  bravais_TYP *G,
              *G_tr,
              *H,
              *H_tr,
              *G_neu,
              *Aut;

  matrix_TYP *erg,
            **base,
            **N,
            **forms,
             *SV,
             *A,
             *id;

  FILE *debug;

  bahn **strong;

  char comment[1000];

  int i,
      tmp;

  read_header(argc, argv);
  if ((FILEANZ < 2) || (is_option('h') && optionnumber('h') ==0)){
     printf("Usage: %s 'file1' 'file2' [-f]\n",argv[0]);
     printf("\n");
     printf("file1: bravais_TYP containing G.\n");
     printf("file1: bravais_TYP containing H.\n");
     printf("\n");
     printf("Searches for a matrix X in GL_n(Z) such that { X g X^-1 | g  in G } = H.\n");
     printf("If no such matrix exists, an error messages is printed to stdout.\n");
     printf("\n");
     printf("WARNING: THE GROUPS IN file1 AND file2 HAVE TO BE FINITE!\n");
     printf("\n");
     printf("Options:\n");
     printf("-f   : recalculate the formspace of the groups even if it\n");
     printf("       is given already.\n");
     printf("-d   : some special files for debugging are created. Do not\n");
     printf("       use this option unless you know what you are doing.\n");
     printf("\n");
     printf("Cf. Is_finite, Bravais_equiv.\n");
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


  /* paranoia setting: recalculate the formspace because it has to be a
     Z-basis */
  if ((G->form !=NULL) && is_option('f')){
     for (i=0;i<G->form_no;i++){
        free_mat(G->form[i]);
     }
     free(G->form);
     G->form = NULL;
  }

  if (G->form == NULL){
    G->form = formspace(G->gen,G->gen_no,1,&tmp);
    G->form_no = tmp;
  }

  if ((H->form !=NULL) && (is_option('f'))){
     for (i=0;i<H->form_no;i++){
        free_mat(H->form[i]);
     }
     free(H->form);
     H->form = NULL;
  }

  if (H->form == NULL){
    H->form = formspace(H->gen,H->gen_no,1,&tmp);
    H->form_no = tmp;
  }

  G_tr = tr_bravais(G,1,FALSE);
  H_tr = tr_bravais(H,1,FALSE);

  if (is_option('d')){ printf("vor is_z_equivalent\n");}

  /* look whether the bravais groups have the same type */
  erg = is_z_equivalent(G,G_tr,H,H_tr);

  if (is_option('d')){
     printf("nach is_z_equivalent\n");
     put_mat(erg,"erg","erg",2);
  }

  if (erg != NULL){
     /* conjugate the groups with this element */
     G_neu = konj_bravais(G,erg);

     /* G_neu is not nessesarily the whole bravais_group */
     /* and apply a trick to speed up for large generating sets */
     id = init_mat(G->dim,G->dim,"1");
     A = rform(G_neu->gen,G_neu->gen_no,id,101);
     free_mat(id);


     forms = (matrix_TYP **) malloc((G->form_no+1)*sizeof(matrix_TYP *));
     forms[0] = A;
     SV = short_vectors(A,max_diagonal_entry(A),0,0,0,&i);
     for (i=0;i<G_neu->form_no;i++){forms[i+1] = G_neu->form[i];}
     Aut = autgrp(forms,G_neu->form_no+1,SV,NULL,0,NULL);
     free_mat(SV);
     free_mat(A);
     free(forms);

     /* stick the right elements into N, i.e. the generators of Aut and
        the generators of the normlizer of G_neu */
     N = (matrix_TYP **) malloc((G_neu->normal_no + Aut->gen_no) *
                                 sizeof(matrix_TYP *));
     for (i=0;i<(G_neu->normal_no+Aut->gen_no);i++){
        if (i<G_neu->normal_no){
           N[i] = G_neu->normal[i];
        }
        else{
           N[i] = Aut->gen[i-G_neu->normal_no];
        }
     }


     /* get strong generators for G_neu */
     base = get_base(G_neu);
     strong = strong_generators(base,G_neu,0);
     G_neu->order = size(strong);
     if (is_option('d')) put_bravais(G_neu,NULL,NULL);

     /* let's see whether we can conjugate G_neu to a be a supergroup
        of H */
     A = conjugated(G_neu,H,N,G_neu->normal_no + Aut->gen_no, strong);

     /* we will need this space for strong generators for H */
     for (i=0;i<G->dim;i++){
        free_bahn(strong[i]);
        free(strong[i]);
     }
     free(strong);

     if (A==NULL){
        /* the groups have the same bravais type, but are not conjugates */
        free_mat(erg);
        erg = NULL;
     }
     else{
        /* they are conjugate iff the have the same order */
        strong = strong_generators(base,H,0);
        if (size(strong) == G_neu->order){
           mat_muleq(A,erg);
           free_mat(erg);
           erg = A;
           A = NULL;
        }
        else{
           free_mat(erg);
           erg = NULL;
           free_mat(A);
        }
        for (i=0;i<G->dim;i++){
           free_bahn(strong[i]);
           free(strong[i]);
        }
        free(strong);
     }

     /* cleaning up */
     free_bravais(G_neu);
     free(N);
     free_bravais(Aut);
     for (i=0;i<G->dim;i++) free_mat(base[i]);
     free(base);
  }

  if (erg == NULL){
     printf("the groups are not conjugated in GL_n(Z)\n");
  }
  else{

     if (is_option('d')){
        debug = fopen("z_equiv.tmp","w");
        fprintf(debug,"#%d\n",G->gen_no + H->gen_no);
        A = long_mat_inv(erg);
        for (i=0;i<G->gen_no;i++){
           SV = mat_mul(erg,G->gen[i]);
           mat_muleq(SV,A);
           fput_mat(debug,SV,NULL,2);
           free_mat(SV);
        }
        free_mat(A);
        for (i=0;i<H->gen_no;i++){
           fput_mat(debug,H->gen[i],NULL,2);
        }
           
        fput_mat(debug,erg,NULL,2);
        fclose(debug);
     }

     sprintf(comment,"conjugates the group of %s to the group of %s",
                      FILENAMES[0],FILENAMES[1]);
     put_mat(erg,NULL,comment,2);
     free_mat(erg);
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
