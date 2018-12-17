#include"typedef.h"
#include"bravais.h"
#include"matrix.h"
#include"symm.h"
#include"autgrp.h"
#include"voronoi.h"

/*******************************************************************************
@
@-------------------------------------------------------------------------------
@  FILE: bravais_tools.h
@
@-------------------------------------------------------------------------------
@
*******************************************************************************/

/*******************************************************************************
@
@------------------------------------------------------------------------------
@
@ bravais_TYP *bravais_group(bravais_TYP *H,int flag)
@
@ Calculates the bravais group of H. the result G of this function
@ will have G->gen and G->form assigned.
@ The program relies on H->form to be correct if it is given.
@ If flag>0, then the function will pair reduce the first form.
@ It is save to use this optiuon, but slower in well behaved examples.
@------------------------------------------------------------------------------
@
******************************************************************************/
bravais_TYP *bravais_group(bravais_TYP *H,int flag)
{
   int i,
       F_no;

   matrix_TYP **F,
               *PF,
               *ID,
               *SV;

   bravais_TYP *G;

   ID = init_mat(H->dim,H->dim,"1");
   PF = rform(H->gen,H->gen_no,ID,101);

   if (!flag){
      SV = short_vectors(PF,max_diagonal_entry(PF),0,0,0,&i);
   }

   /* calculate the formspace if not given */
   if (H->form_no == 0){
      H->form = formspace(H->gen,H->gen_no,1,&H->form_no);
   }

   /* F will hold input for autgrp */
   F_no = H->form_no;
   F = (matrix_TYP **) malloc( (F_no+1) * sizeof(matrix_TYP *));
   for (i=0;i<F_no;i++) F[i+1] = copy_mat(H->form[i]);
   F[0] = PF;

   if (!flag){
      /* call autgrp to get the result */
      G = autgrp(F,F_no+1,SV,NULL,0,NULL);
   }
   else{
      G = pr_aut(F,F_no+1,NULL,0,NULL);
   }

   /* copy the forms to G->form */
   G->form = (matrix_TYP **) malloc(F_no * sizeof(matrix_TYP *));
   G->form_no = F_no;
   for (i=0;i<G->form_no;i++) G->form[i] = F[i+1];

   /* clean up */
   free(F);
   free_mat(PF);
   if (!flag) free_mat(SV);
   free_mat(ID);

   return G;   
}

/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ bravais_TYP *copy_bravais(bravais_TYP *H)
@
@ copies all data form H to the result. Nothing is checked.
@----------------------------------------------------------------------------
@
*****************************************************************************/
bravais_TYP *copy_bravais(bravais_TYP *H)
{

  int i;

  bravais_TYP *G;

  G = init_bravais(H->dim);

  if (H->order != 0){
    G->order = H->order;
    memcpy(G->divisors,H->divisors,100*sizeof(int));
  }
  if (H->gen_no != 0){
    G->gen = (matrix_TYP **) malloc(H->gen_no * sizeof(matrix_TYP *));
    for (i=0;i<H->gen_no;i++) G->gen[i] = copy_mat(H->gen[i]);
    G->gen_no = H->gen_no;
  }
  if (H->cen_no != 0){
    G->cen = (matrix_TYP **) malloc(H->cen_no * sizeof(matrix_TYP *));
    for (i=0;i<H->cen_no;i++) G->cen[i] = copy_mat(H->cen[i]);
    G->cen_no = H->cen_no;
  }
  if (H->normal_no != 0){
    G->normal = (matrix_TYP **) malloc(H->normal_no * sizeof(matrix_TYP *));
    for (i=0;i<H->normal_no;i++) G->normal[i] = copy_mat(H->normal[i]);
    G->normal_no = H->normal_no;
  }
  if (H->form_no != 0){
    G->form = (matrix_TYP **) malloc(H->form_no * sizeof(matrix_TYP *));
    for (i=0;i<H->form_no;i++) G->form[i] = copy_mat(H->form[i]);
    G->form_no = H->form_no;
  }
  if (H->zentr_no != 0){
    G->zentr = (matrix_TYP **) malloc(H->zentr_no * sizeof(matrix_TYP *));
    for (i=0;i<H->zentr_no;i++) G->zentr[i] = copy_mat(H->zentr[i]);
    G->zentr_no = H->zentr_no;
  }

  return G;

}
