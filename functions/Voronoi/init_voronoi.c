#include "typedef.h"
#include "utils.h"
#include "voronoi.h"
#include "matrix.h"
#include "bravais.h"
#include "polyeder.h"
#include "tools.h"
#include "getput.h"
#include "bravais.h"
#include "datei.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: init_voronoi.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ voronoi_TYP *init_voronoi()
@
@ allocates storage for a pointer to voronoi_TYP and sets all intergers
@ and all pointers of it to 0 resp. NULL.
@---------------------------------------------------------------------------
@
\**************************************************************************/
voronoi_TYP *
init_voronoi (void)
{
  voronoi_TYP *V;

  V = (voronoi_TYP *)xmalloc(sizeof(voronoi_TYP));
  V->gram = NULL;
  V->SV_no = 0;
  V->min = 0;
  V->pdet = 0;
  V->prime = 0;
  V->vert_no = 0;
  V->vert = NULL;
  V->pol = NULL;
  V->red_inv = NULL;
  V->T = NULL;
  V->Ti = NULL;
  V->Gtrred = NULL;
  V->SVi = NULL;
  V->stab = NULL;
  V->linstab = NULL;
  V->dir_reps = NULL;

  return(V);
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ void clear_voronoi(V)
@ voronoi_TYP *V;
@
@  frees all pointers that are allocated in the the voronoi_TYP *V
@  and sets them to NULL
@  All integers are set 0.
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
clear_voronoi (voronoi_TYP *V)
{
  int i;

  if(V->gram != NULL)
    free_mat(V->gram);
  V->gram = NULL;
  V->SV_no = 0;
  V->min = 0;
  V->pdet = 0;
  V->prime = 0;
  for(i=0;i<V->vert_no;i++)
    free_wall(&V->vert[i]);
  if(V->vert != NULL)
    free(V->vert);
  V->vert = NULL;
  V->vert_no = 0;
  if(V->pol != NULL)
  {
    free_polyeder(V->pol);
    /* deleted 16/197 tilman
    free(V->pol); */
    V->pol = NULL;
  }
  V->pol = NULL;
  if(V->red_inv != NULL)
    free_mat(V->red_inv);
  V->red_inv = NULL;
  if(V->T != NULL)
    free_mat(V->T);
  V->T = NULL;
  if(V->Ti != NULL)
    free_mat(V->Ti);
  V->Ti = NULL;
  if(V->Gtrred != NULL)
  {
     free_bravais(V->Gtrred);
     /* deleted 16/1/97 tilman
     free(V->Gtrred); */
     V->Gtrred = NULL;
  }
  V->Gtrred = NULL;
  if(V->SVi != NULL)
    free_mat(V->SVi);
  V->SVi = NULL;
  if(V->stab != NULL)
  {
     free_bravais(V->stab);
     /* deleted 16/1/97 tilman
     free(V->stab); */
     V->stab = NULL;
  }
  V->stab = NULL;
  if(V->linstab != NULL)
  {
     free_bravais(V->linstab);
     /* deleted 16/1/97 tilman
     free(V->linstab); */
     V->linstab = NULL;
  }
  V->linstab = NULL;
  if(V->dir_reps != NULL)
    free_mat(V->dir_reps);
  V->dir_reps = NULL;

}



/**************************************************************************\
@---------------------------------------------------------------------------
@ void put_voronoi(V)
@ voronoi_TYP *V;
@
@ Prints the informations in V to standard-output with comments
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
put_voronoi (voronoi_TYP *V)
{
  int i,j;

  put_mat(V->gram, NULL, "matrix of perfect form", 2);
  printf("%d %% number of shortest vectors\n", V->SV_no);
  printf("%d %% minimum of perfect form\n", V->min);
  printf("%d %d %% determinant modulo prime number\n", V->pdet, V->prime);
  printf("%dx%d %% matrix of the Voronoi vertices\n", V->vert_no, V->vert[0]->dim);
  for(i=0;i<V->vert_no;i++)
  {
    for(j=0;j<V->vert[i]->dim;j++)
     printf("%d  ", V->vert[i]->gl[j]);
    printf("\n");
  }
  printf("%% Voronoi polyedral:\n");
  put_polyeder(V->pol);
  put_bravais(V->stab, NULL, "Stabilizer of perfect form");
  put_bravais(V->linstab, NULL, "action of the stabilizer on the space of forms");
  put_mat(V->dir_reps, NULL, "representatives of the Voronoi directions", 0);
}
