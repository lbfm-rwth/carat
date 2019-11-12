#include "typedef.h"
#include "utils.h"
#include "voronoi.h"
#include "longtools.h"
#include "symm.h"
#include "bravais.h"
#include "autgrp.h"
#include "polyeder.h"
#include "matrix.h"
#include "reduction.h"
#include "orbit.h"
#include "tools.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: normalizer.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

static matrix_TYP *
search_normal_isometry (voronoi_TYP *Vneu, voronoi_TYP **V, int Vanz, bravais_TYP *G, bravais_TYP *Gtr, matrix_TYP *bifo, int *no)
{
   int i;
   int found;
   matrix_TYP *X;
 
   found = FALSE;
   *no = -1;
   for(i=0;i<Vanz && found == FALSE;i++)
   {
      X = calc_voronoi_isometry(V[i], Vneu, G, Gtr, bifo);
      if(X != NULL)
      {
       found = TRUE;
       *no = i;
      }
   }
   if(found == TRUE)
   {
     i = *no;
   }
   return(X);
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ voronoi_TYP **normalizer(P, G, Gtr, prime, V_no)
@ matrix_TYP *P;
@ bravais_TYP *G, *Gtr;
@ int prime, *V_no;
@
@  The arguments of normalizer are
@  P:     A G-perfect form (can be calculated with 'first_perfect')
@  G:     A group given as 'bravais_TYP*' where 'G->form' must contain a
@         Z-Basis of the integral G-invariant matrices.
@  Gtr:   The group G^{tr} given as 'bravais_TYP', also
@         with a Z-basis in 'G->form'
@  prime: a prime number (used to calculate determinates modulo p
@         for a simple criterion of two symmetric matrices been
@         not isometric.
@  V_no:  via this pointer the number of returned 'voronoi_TYP'
@         is returned.
@
@ 'normalizer' calculates representatives of the G-perfect forms.
@ The are returned in a list 'voronoi_TYP**'
@ Furthermore generators of the integral normalizer of G are calculated
@ and written to G->normal. Their number is G->normal_no.
@  
@---------------------------------------------------------------------------
@
\**************************************************************************/
voronoi_TYP **
normalizer (matrix_TYP *P, bravais_TYP *G, bravais_TYP *Gtr, int prime, int *V_no)
{
   int i,j,k,l;
   int fdim;
   voronoi_TYP **V;
   voronoi_TYP *Vneu;
   matrix_TYP *Dir, **N, *X, *bifo;
   int Vanz, Nanz, diranz, no;
   int lc, rc;

   fdim = G->form_no;
   V = (voronoi_TYP **)xmalloc(1 *sizeof(voronoi_TYP *));
   V[0] = init_voronoi();
   V[0]->gram = copy_mat(P);
   Check_mat(V[0]->gram);
   bifo = trace_bifo(G->form, Gtr->form, fdim);
   calc_voronoi_complete(V[0], G, Gtr, bifo, prime);
   Nanz = V[0]->stab->gen_no;
   N = (matrix_TYP **)xmalloc(Nanz *sizeof(matrix_TYP *));
   for(i=0;i<Nanz;i++)
     N[i] = copy_mat(V[0]->stab->gen[i]);
   Vanz = 1;
   Vneu = init_voronoi();
   for(i=0;i<Vanz;i++)
   {
      diranz = V[i]->dir_reps->cols;
      real_mat(V[i]->dir_reps, 3, diranz);

      /* geaendert testweise am 06.11.96 tilman */
      for(j=0;j<diranz;j++)
      {
         /***************************************************************\
         | Calculate 'Vneu', the neighbour of V[i] in direction j
         \***************************************************************/
         Dir = vec_to_form(V[i]->pol->vert[V[i]->dir_reps->array.SZ[0][j]]->v,
                           G->form, fdim);
         Vneu->gram = voronoi_neighbour(V[i]->gram, Dir, V[i]->min, &lc, &rc);
         free_mat(Dir);

         /* changed: tilman 07/11/96 */
         /* voronoi_neighbour returns NULL sometimes, and this causes
            chrashes in calc_voronoi_basics */
         if (Vneu->gram != NULL){
            divide_by_gcd(Vneu->gram);
            calc_voronoi_basics(Vneu, G, Gtr, prime);
            /***************************************************************\
            | Check if 'Vneu' is a new typ of perfect form or not
            | if 'Vneu' is new, 'no = -1' and 'X = NULL'
            | if not, 'X' is an isometry in the integral normalizer of 'G'
            |          and 'Vneu' is isometric to V[no].
            \***************************************************************/
            X = search_normal_isometry(Vneu, V, Vanz, G, Gtr, bifo, &no);
            if(no == -1)
            {
              V[i]->dir_reps->array.SZ[2][j] = Vanz;
              calc_voronoi_complete(Vneu, G, Gtr, bifo, prime);

              l = Vneu->stab->gen_no;
              Nanz += l;
              N =(matrix_TYP **)xrealloc(N,Nanz *sizeof(matrix_TYP *));
              for(k=0;k<l; k++)
                N[Nanz-l+k] = copy_mat(Vneu->stab->gen[k]);

              Vanz++;
              V=(voronoi_TYP **)xrealloc(V,Vanz *sizeof(voronoi_TYP *));
              V[Vanz-1] = Vneu;
              Vneu = init_voronoi();
            }
            else
            {
               V[i]->dir_reps->array.SZ[2][j] = no;
               Nanz ++;
               N=(matrix_TYP **)xrealloc(N,Nanz *sizeof(matrix_TYP *));
               N[Nanz-1] =  X;
               clear_voronoi(Vneu);
            }
         }
      }
   }
   *V_no = Vanz;
   if(G->normal_no != 0)
   {
     for(i=0;i<G->normal_no;i++)
       free_mat(G->normal[i]);
     free(G->normal);
   }
   G->normal = N;
   G->normal_no = Nanz;
   free_mat(bifo);
   if (Vneu != NULL){
      clear_voronoi(Vneu);
      free(Vneu);
      Vneu = NULL;
   }

   /* inserted tilman 20.06.97 to reduce the number of generators for
   the normalizer */
   red_normal(G);

   return(V);
}
