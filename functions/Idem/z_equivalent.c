#include "typedef.h"
#include "getput.h"
#include "name.h"
#include "gmp.h"
#include "bravais.h"
#include "datei.h"
#include "matrix.h"
#include "voronoi.h"
#include "autgrp.h"
#include "symm.h"
#include "base.h"
#include "zass.h"
#include "longtools.h"
#include "idem.h"

matrix_TYP *z_equivalent(bravais_TYP *G,
                         bravais_TYP **G_tr,
                         bravais_TYP *H)
{

  bravais_TYP *H_tr,
              *G_neu,
              *Aut;

  matrix_TYP *erg,
            **base,
            **N,
            **forms,
             *SV,
             *A,
             *id;

  bahn **strong;

  int i,
      tmp;

  if (G->form == NULL){
    G->form = formspace(G->gen,G->gen_no,1,&tmp);
    G->form_no = tmp;
  }

  if (H->form == NULL){
    H->form = formspace(H->gen,H->gen_no,1,&tmp);
    H->form_no = tmp;
  }

  if (*G_tr == NULL)
      *G_tr = tr_bravais(G,1,FALSE);

  H_tr = tr_bravais(H,1,FALSE);


  /* look whether the bravais groups have the same type */
  erg = is_z_equivalent(G,*G_tr,H,H_tr);

  if (INFO_LEVEL & 4){
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
     strong = strong_generators(base,G_neu,FALSE);
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
        strong = strong_generators(base,H,FALSE);
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

  free_bravais(H_tr);

  return erg;

}

