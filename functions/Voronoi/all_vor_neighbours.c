#include "typedef.h"
#include "utils.h"
#include "matrix.h"
#include "voronoi.h"
#include "polyeder.h"
#include "tools.h"
#include "bravais.h"

/**********************************************************************\
@
@--------------------------------------------------------------------
@  FILE: all_vor_neighbours.c
@--------------------------------------------------------------------
@
\**********************************************************************/



/**********************************************************************\
@ matrix_TYP *all_voronoi_neighbours(P, G, Ftr, tr_bifo)
@ matrix_TYP *P, **Ftr, *tr_bifo;
@ bravais_TYP *G;
@--------------------------------------------------------------------
@
@ The arguments of 'all_voronoi_neighbours' are:
@ P:		a G-perfect form
@ G:		the group for which P is G-perfect
@ Ftr:		a basis for the integral forms in the formspace of G^{tr}
@ trbifo:	The matrix of the bilinear form
@		F(G)xF(G^{tr})--> R: (A,B)---> trace(AB)
@		with respect to the Z-bases given by G->form and Ftr
@ Let Fdim be the dimension of the formspace of G and m be the number of
@ G-perfect forms that are neighbours of P.
@ Then the result of 'all_perfect_neighbours' is a matrix N with m rows
@ and Fdim+3 columns.
@ For 0 < i< m let R[i] be the form in F(G) represented by the
@ first Fdim entries of the i-th row of N, t.m.
@     R[i] = N[i][0] * G->form[0] + ..... + N[i][Fdim-1] * G->form[Fdim-1]
@ Further let
@     p[i] := N[i][Fdim], r[i] := N[i][Fdim+1], c[i] := N[i][Fdim+2]
@ Then P[i] := (p[i]*P + r[i]*R[i]) / c[i]
@ is a G-perfect form  neighboured to P.
@--------------------------------------------------------------------
@
\**********************************************************************/
matrix_TYP *
all_voronoi_neighbours (matrix_TYP *P, bravais_TYP *G, matrix_TYP **Ftr, matrix_TYP *tr_bifo)
{
  int i,j,k,l, Fdim, Gdim;
  int Pmin, Vanz, Wanz;
  matrix_TYP *erg;
  matrix_TYP **V, *N, *X;
  polyeder_TYP *pol;
  wall_TYP **W;
  int *ww, *pw;
  int pos, lc, rc, g;


  Fdim = G->form_no;
  Gdim = G->form[0]->cols;
  /**********************************************************************\
  | calculate the vertices of the Voronoidomain in F(G^{tr})
  | and the forms in F(G) corresponding to the walls of the Voronoidomain.
  | The coordinates of these forms with respect to the basis given in G->form
  | are obtained as the coordinates of the vetrices of the polyeder_TYP *pol.
  \**********************************************************************/
  V = voronoi_vertices(P, G, &Vanz, &Pmin, &i);
  Wanz = Vanz;
  W = (wall_TYP **)xmalloc(Wanz *sizeof(wall_TYP *));
  ww = (int *)xmalloc(Fdim *sizeof(int));
  pw = (int *)xmalloc(Fdim *sizeof(int));
  for(i=0;i<Vanz;i++)
  {
    W[i] = init_wall(Fdim);
    form_to_vec(ww, V[i], Ftr, Fdim, &k);
    for(j=0;j<Fdim; j++)
      for(k=0;k<Fdim; k++)
         W[i]->gl[j] += ww[k] * tr_bifo->array.SZ[j][k];
    normal_wall(W[i]);
    free_mat(V[i]);
  }
  free(V);
  pol = first_polyeder(W, Wanz);
  if(pol == NULL)
  {
     printf("Error in all_perfect_neighbours:  P not G-perfekt\n");
     exit(3);
  }
  for(i=0;i<Wanz;i++)
  {
    refine_polyeder(pol, W[i]);
    free_wall(&W[i]);
  }
  free(W);

  /*
  put_polyeder(pol);
  */

  /******************************************************************\
  | Each G-perfect form A, that is a neighbour of P, is (up to
  | multiplication with a positive scalar) given by
  | N = lc *P + rc * X,
  | where X is a form defined by a vertex of the polyeder pol.
  | The coefficients are calculated with the function 'voronoi_neighbour'. 
  \******************************************************************/
  form_to_vec(pw, P, G->form, Fdim, &k);
  erg = init_mat(pol->vert_no, (Fdim+3), "");
  X = init_mat(Gdim, Gdim, "");
  X->flags.Integral = X->flags.Symmetric = TRUE;
  X->flags.Diagonal = FALSE;
  pos = 0;
  for(i=0;i<pol->vert_no;i++)
  {
    /******************************************************\
    | Calculate X
    \******************************************************/
    for(j=0;j<Gdim;j++)
      for(k=0;k<=j;k++)
      {
        X->array.SZ[j][k] = 0;
        for(l=0;l<Fdim;l++)
          X->array.SZ[j][k] +=pol->vert[i]->v[l] * G->form[l]->array.SZ[j][k];
        X->array.SZ[k][j] = X->array.SZ[j][k];
      }
    N = voronoi_neighbour(P, X, Pmin, &lc, &rc);
    /******************************************************\
    | Divide N by the greatest common divisor of its entries
    \******************************************************/
    if(N != NULL)
    {
       for(j=0;j<Fdim;j++)
         erg->array.SZ[pos][j] = pol->vert[i]->v[j];
       erg->array.SZ[pos][Fdim] = lc;
       erg->array.SZ[pos][Fdim+1] = rc;
       for(j=0;j<Fdim;j++)
         ww[j] = lc * pw[j] + rc * erg->array.SZ[pos][j];
       k = 0;
       while(k<Fdim && ww[k] == 0)
           k++;
       g = ww[k];
       for(j=k+1;j<Fdim;j++)
       {
         if(ww[j] != 0)
         {
           l = GGT(g, ww[j]);
           g = l;
         }
       }
       erg->array.SZ[pos][Fdim+2] = g;
       free_mat(N);
       pos++;
    }
  }

  /* changed on 14/1/97 tilman from
  erg->rows = pos;
  to: */
  real_mat(erg,pos,erg->cols);

  free_mat(X);
  free(pw); free(ww);
/*********************************
  put_polyeder(pol);
***********************************/
  free_polyeder(pol);
  return(erg);
}
