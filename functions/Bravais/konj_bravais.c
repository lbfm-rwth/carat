#include"typedef.h"
#include"matrix.h"
#include"longtools.h"
#include"bravais.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: konj_bravais.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ bravais_TYP *konj_bravais(B, T)
@ bravais_TYP *B;
@ matrix_TYP *T;
@
@  calculates the group T B T^(-1)
@  all informations of B are transformed: B->gen, B->form,....
@---------------------------------------------------------------------------
@
\**************************************************************************/


bravais_TYP *
konj_bravais (bravais_TYP *B, matrix_TYP *T)
{
   int i,dim;
   bravais_TYP *G;
   matrix_TYP *Ti, *Titr, *waste;

   G = init_bravais(B->dim);
   Ti = long_mat_inv(T);

   /* inserted tilman 11/4/97 */
   if (B->order != 0){
     G->order = B->order;
     memcpy(G->divisors,B->divisors,100*sizeof(int));
   }

   if(B->gen_no != 0)
   {
     G->gen_no = B->gen_no;
     if( (G->gen = (matrix_TYP **)malloc(B->gen_no *sizeof(matrix_TYP *))) == NULL)
     {
        printf("malloc of 'G->gen' in 'konj_bravais' failed\n");
        exit(2);
     }
     for(i=0;i<B->gen_no;i++)
     {
        waste = mat_mul(T, B->gen[i]);
        G->gen[i] = mat_mul(waste, Ti);
        free_mat(waste);
     }
   }

   if(B->form_no != 0)
   {
     Titr = tr_pose(Ti);
     G->form_no = B->form_no;
     if( (G->form = (matrix_TYP **)malloc(B->form_no *sizeof(matrix_TYP *))) == NULL)
     {
        printf("malloc of 'G->form' in 'konj_bravais' failed\n");
        exit(2);
     }
     for(i=0;i<B->form_no;i++)
     {
        waste = mat_mul(Titr, B->form[i]);
        G->form[i] = mat_mul(waste, Ti);
        free_mat(waste);
     }
     free_mat(Titr);
     /* don't do this, it causes trouble in normalizer because
     normalizer assumes the space of forms is calculated in the above way:
     long_rein_formspace(G->form,G->form_no,1); */
   }

   if(B->zentr_no != 0)
   {
     G->zentr_no = B->zentr_no;
     if( (G->zentr = (matrix_TYP **)malloc(B->zentr_no *sizeof(matrix_TYP *))) == NULL)
     {
        printf("malloc of 'G->zentr' in 'konj_bravais' failed\n");
        exit(2);
     }
     for(i=0;i<B->zentr_no;i++)
     {
        waste = mat_mul(T, B->zentr[i]);
        G->zentr[i] = mat_mul(waste, Ti);
        free_mat(waste);
     }
   }

   if(B->normal_no != 0)
   {
     G->normal_no = B->normal_no;
     if( (G->normal = (matrix_TYP **)malloc(B->normal_no *sizeof(matrix_TYP *))) == NULL)
     {
        printf("malloc of 'G->normal' in 'konj_bravais' failed\n");
        exit(2);
     }
     for(i=0;i<B->normal_no;i++)
     {
        waste = mat_mul(T, B->normal[i]);
        G->normal[i] = mat_mul(waste, Ti);
        free_mat(waste);
     }
   }

   if(B->cen_no != 0)
   {
     G->cen_no = B->cen_no;
     if( (G->cen = (matrix_TYP **)malloc(B->cen_no *sizeof(matrix_TYP *))) == NULL)
     {
        printf("malloc of 'G->cen' in 'konj_bravais' failed\n");
        exit(2);
     }
     for(i=0;i<B->cen_no;i++)
     {
        waste = mat_mul(T, B->cen[i]);
        G->cen[i] = mat_mul(waste, Ti);
        free_mat(waste);
     }
   }
  
   free_mat(Ti);
   return(G);
}
