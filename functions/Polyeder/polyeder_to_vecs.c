#include "typedef.h"
#include "utils.h"
#include"matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: polyeder_to_vecs.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP **polyeder_to_vecs(P)
@ polyeder_TYP *P;
@
@ polyeder_to_vecs calculates two matrices X[0] and X[1]
@ where the rows of X[0] contain the coordinates of the vertices
@ i.e. X[0]->array.SZ[i] = P->vert[i]->v
@ and the rows of X[1] contain the coordinates of the walls
@ i.e. X[1]->array.SZ[i] = P->wall[i]->gl.
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP **
polyeder_to_vecs (polyeder_TYP *P)
{
  int i,j, dim;
  matrix_TYP **M;

  M = (matrix_TYP **)xmalloc(2 *sizeof(matrix_TYP *));
  dim = P->vert[0]->dim;
  M[0] = init_mat(P->vert_no, dim, "");
  M[1] = init_mat(P->wall_no, dim, "");
  for(i=0;i<P->vert_no;i++)
    for(j=0;j<dim;j++)
      M[0]->array.SZ[i][j] = P->vert[i]->v[j];
  for(i=0;i<P->wall_no;i++)
    for(j=0;j<dim;j++)
      M[1]->array.SZ[i][j] = P->wall[i]->gl[j];
  return(M);
}
