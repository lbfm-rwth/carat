#include "typedef.h"
#include "utils.h"
#include"matrix.h"
#include"bravais.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  tr_bravais.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ bravais_TYP *tr_bravais(B, calcforms, invert)
@ bravais_TYP *B;
@ int calcforms;
@ int invert;
@
@  'tr_bravais' caclulates the group G =  B^{tr}
@  The matrices in B->gen, B->zentr, B->normal and B->cen are transposed
@  (if exists)
@  If calcforms == 0, then the matrices of G->form are not calculated
@  If calcforms == 1, then G->form is calculated with 'formspace'
@  If calcforms == 2, then G->form is calculated with 'invar_space'
@                     if B->form_no != 0, invar_space is started with
@                     fdim = B->form_no, otherwise the argument 'fdim'
@                     for invar_space is calculated by calculating the
@                     formspace with p_formspace modulo 101.
@  If (invert == TRUE) the matrices in G->gen,G->cen,G->normal
@  are inverted.
@---------------------------------------------------------------------------
@
\**************************************************************************/
bravais_TYP *
tr_bravais (bravais_TYP *B, int calcforms, int invert)
{
   bravais_TYP *G;
   int i,j;
   matrix_TYP **XX,
               *tmp;

   G = init_bravais(B->dim);
   if(B->gen_no != 0)
   {
     G->gen_no = B->gen_no;
     G->gen = (matrix_TYP **)xmalloc(G->gen_no *sizeof(matrix_TYP *));
     for(i=0;i<G->gen_no;i++)
       if (invert){
          tmp = mat_inv(B->gen[i]);
          G->gen[i] = tr_pose(tmp);
          free_mat(tmp);
       }
       else{
          G->gen[i] = tr_pose(B->gen[i]);
       }
   }
   if(B->zentr_no != 0)
   {
     G->zentr_no = B->zentr_no;
     G->zentr = (matrix_TYP **)xmalloc(G->zentr_no *sizeof(matrix_TYP *));
     for(i=0;i<G->zentr_no;i++)
       G->zentr[i] = tr_pose(B->zentr[i]);
   }
   if(B->normal_no != 0)
   {
     G->normal_no = B->normal_no;
     G->normal = (matrix_TYP **)xmalloc(G->normal_no *sizeof(matrix_TYP *));
     for(i=0;i<G->normal_no;i++)
       if (invert){
          tmp = mat_inv(B->normal[i]);
          G->normal[i] = tr_pose(tmp);
          free_mat(tmp);
       }
       else{
          G->normal[i] = tr_pose(B->normal[i]);
       }
   }
   if(B->cen_no != 0)
   {
     G->cen_no = B->cen_no;
     G->cen = (matrix_TYP **)xmalloc(G->cen_no *sizeof(matrix_TYP *));
     for(i=0;i<G->cen_no;i++)
       if (invert){
          tmp = mat_inv(B->cen[i]);
          G->cen[i] = tr_pose(tmp);
          free_mat(tmp);
       }
       else{
          G->cen[i] = tr_pose(B->cen[i]);
       }
   }
   if(calcforms == 1)
       G->form = formspace(G->gen, G->gen_no, 1, &G->form_no);
   if(calcforms == 2)
   {
       if(B->form_no != 0)
          G->form = invar_space(G->gen, G->gen_no, B->form_no, 1, 100, &G->form_no); 
       else
       {
          XX = p_formspace(G->gen, G->gen_no, 101, 1, &i);
          for(j=0;j<i;j++)
            free_mat(XX[i]);
          free(XX);
          G->form = invar_space(G->gen, G->gen_no, i, 1, 100, &G->form_no); 
       }
   }

   /* inserted a return: tilman 7/1/97 */
   return G;
}
