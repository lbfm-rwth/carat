#include "typedef.h"
#include "matrix.h"
#include "longtools.h"
#include "sort.h"
#include "bravais.h"

static matrix_TYP *formspace_operation(matrix_TYP **F,int Fno,matrix_TYP *N)
{

   int i;

   matrix_TYP *res,
              *tmp,
              *N_tr;

   res = init_mat(Fno,Fno,"");
   N_tr = tr_pose(N);

   for (i=0;i<Fno;i++){
      tmp = scal_pr(N_tr,F[i],TRUE);
      form_to_vec_modular(res->array.SZ[i],tmp,F,Fno);
      free_mat(tmp);
   }

   free_mat(N_tr);

   return res;
}


static int position(matrix_TYP **a,matrix_TYP *x,int n)
/* returns the first index i<n such that a[i] == x, and
   -1 if none exists */
{
  int i=0;

  while (i<n){
    if (mat_comp(a[i],x) == 0){
       return i;
    }
    i++;
  }
  return -1;
}

void red_normal(bravais_TYP *G)
{

   int i;

   matrix_TYP **REP,
               *tmp;

  REP = (matrix_TYP **) malloc(G->normal_no * sizeof(matrix_TYP *));

   /* calculate the presentation on the formspace (bare in mind that it
      is faithfull for N_GL_n(Z) (G)/G if G is a bravais_group */
   for (i=0;i<G->normal_no;i++){
      REP[i] = formspace_operation(G->form,G->form_no,G->normal[i]);
   }

   /* see if the are nessesary */
   for (i=1;i<G->normal_no;i++){
      if (position(REP,REP[i],i) != (-1)){
         /* throw it away */
         free_mat(G->normal[i]);
         G->normal[i] = NULL;
      }
      else{
         /* we might got an inverse already */
         tmp = long_mat_inv(REP[i]);
         if (position(REP,tmp,i) != (-1)){
            free_mat(G->normal[i]);
            G->normal[i] = NULL;
         }
         free_mat(tmp);
      }
   }

   /* free REP */
   for (i=0;i<G->normal_no;i++)
      free_mat(REP[i]);
   free(REP);

   /* now swap out the NULL's we got in G->normal */
   for (i=0;i<G->normal_no;i++){
      if (G->normal[i] == NULL){
         G->normal_no--;
         G->normal[i]= G->normal[G->normal_no];
         i--;
      }
   }

   return;
}
