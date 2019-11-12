#include "typedef.h"
#include "utils.h"
#include"voronoi.h"
#include"longtools.h"
#include"symm.h"
#include"bravais.h"
#include"autgrp.h"
#include"polyeder.h"
#include"matrix.h"
#include"reduction.h"
#include"orbit.h"
#include"datei.h"
#include"tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: calc_voronoi_data.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ void calc_voronoi_basics(V, G, Gtr, prime)
@ voronoi_TYP *V;
@ bravais_TYP *G, *Gtr;
@ int prime;
@
@ For using this function V->gram must be a G-perfect form
@ The function calculates:
@
@  V->SV_no:       the number of shortest vectors of V->gram.
@  V->min:         the minimum of V->gram.
@  V->pdet:        the determinante of V->gram modulo prime.
@  V->prime:       = prime.
@  V->vert:        the coordinates of the voronoi vertices with
@                  respect to the Z-basis of the space of invariant
@                  symmetric matrices of G^{tr} given in Gtr->form.
@  V->vert_no:     the number of voronoi_vertices.
@  V->red_inv:     the pair-reduced of (V->gram)^{-1} (with kgv = 1).
@  V->T:           the matrix such that
@                     T  *(V->gram)^{-1} * T^{tr} = V->red_inv
@  V->Ti:          the inverse matrix of V->T
@  V->Gtrred:      the bravais_TYP obtained by konjugating Gtr with
@                  the transposed of V->Ti with the function
@                  'konj_bravais'
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
calc_voronoi_basics (voronoi_TYP *V, bravais_TYP *G, bravais_TYP *Gtr, int prime)
{
  int i,k;
  int anz;
  matrix_TYP **XX, *TTri;

  XX = voronoi_vertices(V->gram, G, &anz, &V->min, &V->SV_no);
  if(V->vert_no != 0 && V->vert != NULL)
  {
     for(i=0;i<V->vert_no;i++)
       free_wall(&V->vert[i]);
     free(V->vert);
  }
  V->vert_no = anz;
  V->vert = (wall_TYP **)xmalloc(anz *sizeof(wall_TYP *));
  for(i=0; i<anz;i++)
  {
    V->vert[i] = init_wall(G->form_no);
    form_to_vec(V->vert[i]->gl, XX[i], Gtr->form, Gtr->form_no, &k);
    free_mat(XX[i]);
  }
  free(XX);

  V->pdet = p_mat_det(V->gram, prime);

  V->prime = prime;
  if(V->red_inv == NULL)
  {
    V->T = init_mat(V->gram->cols, V->gram->cols, "");
    for(i=0;i<V->gram->cols;i++)
       V->T->array.SZ[i][i] = 1;
    V->red_inv = pair_red_inv(V->gram, V->T);
    V->Ti = long_mat_inv(V->T);
    TTri = tr_pose(V->Ti);
    V->Gtrred = konj_bravais(Gtr, TTri);
    free_mat(TTri);
  }
}




/**************************************************************************\
@---------------------------------------------------------------------------
@ void calc_voronoi_pol(V, bifo)
@ voronoi_TYP *V;
@
@  calculates the polyeder_TYP V->pol where the walls of V->pol
@  are the walls in 'V->vert' multiplied with the transposed of bifo,
@  i.e. V->vert[i] * bifo^{tr} 0 <= i <V->vert_no
@  is the input for 'first_polyeder' and 'refine_polyeder' used
@ to calculate 'V->pol'.
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
calc_voronoi_pol (voronoi_TYP *V, matrix_TYP *bifo)
{
   int i,j,k;
   int fdim;
   wall_TYP **W;

   fdim = bifo->cols;
  
   W = (wall_TYP **)xmalloc(V->vert_no *sizeof(wall_TYP *));
   for(i=0;i<V->vert_no;i++)
   {
       W[i] = init_wall(fdim);
       for(j=0;j<fdim;j++)
        for(k=0;k<fdim;k++)
         W[i]->gl[j] += V->vert[i]->gl[k] * bifo->array.SZ[j][k];
   }
   V->pol = first_polyeder(W, V->vert_no);
   for(i=0;i<V->vert_no;i++)
   {
     j = refine_polyeder(V->pol, W[i]);

     /* tilman : changed 17/1/97 from 
     if(j == FALSE){
      free_wall(&W[i]);
     }  to: */
     free_wall(&W[i]);
   }
   free(W);
}




/**************************************************************************\
@---------------------------------------------------------------------------
@ void calc_voronoi_good_inv(V, Gtr)
@ voronoi_TYP *V;
@ bravais_TYP *Gtr;
@
@ applies 'short_red' and 'mink_red' to V->red_inv to obtain a better
@ form.
@ V->T, V->Ti and V->Gtrred are changed in the same manner as described
@ for 'calc_voronoi_basics'
@
@ Furthermore the shortest vectors of the better reduced form are
@ calculated up to the norm that equals the maximal diagonal entry
@ of V->red_inv.
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
calc_voronoi_good_inv (voronoi_TYP *V, bravais_TYP *Gtr)
{
   int i,j,maxd;
   int dim, mr;
   matrix_TYP *A, *T, *T1;
   int has_changed;

   dim = V->gram->cols;
   has_changed = FALSE;
  
   T = init_mat(dim, dim, "");
   for(i=0;i<dim;i++)
     T->array.SZ[i][i] = 1;
   A = pr_short_red(V->red_inv, T);
   maxd = V->red_inv->array.SZ[0][0];
   for(i=1;i<dim;i++)
   { if(V->red_inv->array.SZ[i][i] > maxd)
       maxd = V->red_inv->array.SZ[i][i];
   }
   j = A->array.SZ[0][0];
   for(i=1;i<dim;i++)
   { if(A->array.SZ[i][i] > j)
       j = A->array.SZ[i][i];
   }
   if(j<=maxd)
   {
     has_changed = TRUE;
     T1 = mat_mul(T, V->T);
     free_mat(V->T);
     V->T = T1;
     T1 = NULL;
     free_mat(T);
     T = NULL;
     maxd = j;
     free_mat(V->red_inv);
     V->red_inv = A;
     A = NULL;
   }
   else
     free_mat(A);
   mr = TRUE;
   for(i=0;i<dim;i++)
   {
      if(V->red_inv->array.SZ[i][i] < maxd)
         mr = FALSE;
   }
   if(mr == FALSE)
   {
     if(T == NULL)
       T = init_mat(dim, dim, "1");
     for(i=0;i<dim;i++)
     {
       for(j=0;j<dim;j++)
         T->array.SZ[i][j] = 0;
       T->array.SZ[i][i] = 1;
     }
     /* put_mat(V->red_inv, NULL, "veredinv vor minkred", 2); */
     A = mink_red(V->red_inv, T);
     /* put_mat(V->red_inv, NULL, "veredinv nach minkred", 2);
     put_mat(A, NULL, "A nach minkred", 2); */
     j = A->array.SZ[0][0];
     for(i=1;i<dim;i++)
     { if(A->array.SZ[i][i] > j)
         j = A->array.SZ[i][i];
     }
     if(j<=maxd)
     {
       has_changed = TRUE;
       T1 = mat_mul(T, V->T);
       free_mat(V->T);
       V->T = T1;
       T1 = NULL;
       maxd = j;
       free_mat(V->red_inv);
       V->red_inv = A;
       A = NULL;
     }
     else
       free_mat(A);
   }

   /* changed on 14/1/97 tilman: (inserted if (..)) */
   if (T!=NULL){
      free_mat(T);
   }

   if(has_changed == TRUE)
   {
     free_mat(V->Ti);
     V->Ti = long_mat_inv(V->T);
     T = tr_pose(V->Ti);
     free_bravais(V->Gtrred);
     V->Gtrred = konj_bravais(Gtr, T);
     free_mat(T);
   }
   maxd = V->red_inv->array.SZ[0][0];
   for(i=0;i<dim;i++)
   {
     if(V->red_inv->array.SZ[i][i] > maxd)
        maxd = V->red_inv->array.SZ[i][i];
   }
   V->SVi = short_vectors(V->red_inv, maxd, 0, 0, 0, &i);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ void calc_voronoi_stab(V, G, Gtr, bifo)
@ voronoi_TYP *V;
@ bravais_TYP *G, *Gtr;
@ matrix_TYP *bifo;
@
@ calculates a generating set for
@  H = { h\in GL_n(Z) | h^{tr} * V->gram * h = V->gram and
@                       h G h^{-1} = G                      }
@
@  Furthermore the bravais_TYP* V->linstab is calculated, that
@  describes the linear action for H on the space of G-invariant
@  symmetric matrices with respect to the Z-basis given in 
@  G->form (c.f. the function 'nomlin' in directory 'Bravais'
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
calc_voronoi_stab (voronoi_TYP *V, bravais_TYP *G, bravais_TYP *Gtr, matrix_TYP *bifo)
{

   int i,j;
   int dim,fdim, Vanz;
   matrix_TYP *M, **bas, *M1, **W;
   bravais_TYP *S;

   dim = V->gram->cols;
   fdim = bifo->cols;

   if(V->prime == 0)
     calc_voronoi_basics(V, G, Gtr, 101);
   Vanz = V->vert_no;
   if(V->SVi == NULL)
     calc_voronoi_good_inv(V, Gtr);
   /*******************************************************************\
   | calculate 'bas', a consisting of vertices given as symetric matrices
   \*******************************************************************/
   M = init_mat(Vanz, fdim, "");
   for(i=0;i<Vanz;i++)
     for(j=0;j<fdim;j++)
       M->array.SZ[i][j] = V->vert[i]->gl[j];
   M1 = long_qbase(M);
   free_mat(M);
   if(M1->rows != fdim)
   {
     printf("error in 'calc_voronoi_stab': form not G-perfect\n");
     exit(3);
   }
   bas = (matrix_TYP **)xmalloc(fdim *sizeof(matrix_TYP *));
   for(i=0;i<fdim;i++)
     bas[i] = vec_to_form(M1->array.SZ[i], V->Gtrred->form, fdim);
   free_mat(M1);

   /*******************************************************************\
   | calculate 'W', the symmetric matrices given by the vertices in V->vert
   \*******************************************************************/
   W = (matrix_TYP **) xmalloc(Vanz *sizeof(matrix_TYP *));
   for(i=0;i<Vanz;i++)
     W[i] = vec_to_form(V->vert[i]->gl, V->Gtrred->form, fdim);

   /*******************************************************************\
   | calculate  the stabilizer
   \*******************************************************************/
   S = perfect_normal_autgrp(V->red_inv, V->SVi, V->Gtrred->gen, Gtr->gen_no, NULL, W, Vanz, bas, fdim);
   for(i=0;i<fdim;i++)
     free_mat(bas[i]);
   free(bas);
   for(i=0;i<Vanz;i++)
     free_mat(W[i]);
   free(W);

   V->stab = init_bravais(dim);
   V->stab->gen_no = S->gen_no;
   V->stab->gen = (matrix_TYP **)xmalloc(S->gen_no *sizeof(matrix_TYP *));
   for(i=0;i<S->gen_no;i++)
   {
     M = tr_pose(S->gen[i]);
     M1 = mat_mul(V->Ti, M);
     V->stab->gen[i] = mat_mul(M1, V->T);
     free_mat(M); free_mat(M1);
   }
   V->stab->order = S->order;
   for(i=0;i<100;i++)
    V->stab->divisors[i] = S->divisors[i];
   free_bravais(S);
   V->linstab = init_bravais(fdim);
   V->linstab->gen = (matrix_TYP **)xmalloc(V->stab->gen_no *sizeof(matrix_TYP *));
   V->linstab->gen_no = V->stab->gen_no;
   for(i=0;i<V->stab->gen_no;i++)
     V->linstab->gen[i] = normlin(G->form, V->stab->gen[i], fdim);
}




/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *calc_voronoi_isometry(V1, V2, G, Gtr, bifo)
@ voronoi_TYP *V1, *V2;
@ bravais_TYP *G, *Gtr;
@ matrix_TYP *bifo;
@
@  calculates a matrix X in GL_N(Z) such that
@     X * V1->gram * X^{tr} = V2->gram and  X * G * X^{-1} = G
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
calc_voronoi_isometry (voronoi_TYP *V1, voronoi_TYP *V2, bravais_TYP *G, bravais_TYP *Gtr, matrix_TYP *bifo)
{

   int i,j;
   int dim,fdim, Vanz;
   matrix_TYP *M, **bas, *M1, **W;
   matrix_TYP *SV1, *SV2;
   matrix_TYP *X, *E;
   voronoi_TYP *Vn1, *Vn2;
   int md1, md2, md, permute;

   dim = V1->gram->cols;
   fdim = bifo->cols;

   if(V1->prime == 0)
     calc_voronoi_basics(V1, G, Gtr, 101);
   if(V2->prime == 0)
     calc_voronoi_basics(V2, G, Gtr, 101);
   /****************************************************************\
   | make first easy tests, if an isiometry can exist
   \****************************************************************/
    if(V1->prime == V2->prime && V1->pdet != V2->pdet)
      return(NULL);
    if(V1->min != V2->min)
      return(NULL);
    if(V1->SV_no != V2->SV_no)
      return(NULL);
    if(V1->vert_no != V2->vert_no)
      return(NULL);

   Vanz = V1->vert_no;
   if(V1->SVi == NULL)
     calc_voronoi_good_inv(V1, Gtr);

   /*******************************************************************\
   | Check which 'red_inv' has the smaller maximal diagonal entry
   \*******************************************************************/
    md1 = V1->red_inv->array.SZ[0][0];
    md2 = V2->red_inv->array.SZ[0][0];
    for(i=1;i<dim;i++)
    {
       if(V1->red_inv->array.SZ[i][i] > md1)
         md1 = V1->red_inv->array.SZ[i][i];
       if(V2->red_inv->array.SZ[i][i] > md2)
         md2 = V2->red_inv->array.SZ[i][i];
    }
    if(md2 < md1)
    { permute = TRUE; md = md2; Vn1 = V2; Vn2 = V1;}
    else
    { permute = FALSE; md = md1; Vn1 = V1; Vn2 = V2;}

   /*******************************************************************\
   | calculate  the shortest_vectors
   \*******************************************************************/
    if(Vn1->SVi != NULL)
      SV1 = Vn1->SVi;
    else
      SV1 = short_vectors(Vn1->red_inv, md, 0, 0, 0, &i);
    SV2 = short_vectors(Vn2->red_inv, md, 0, 0, 0, &j);
    if(SV1->rows != SV2->rows)
    {
       if(Vn1->SVi == NULL)
         free_mat(SV1);
       free_mat(SV2);
       /* printf("HIER\n"); */
       return(NULL);
    }


   /*******************************************************************\
   | calculate 'bas', a basis consisting of vertices of Vn2
   |  given as symmetric matrices
   \*******************************************************************/
   M = init_mat(Vanz, fdim, "");
   for(i=0;i<Vanz;i++)
     for(j=0;j<fdim;j++)
       M->array.SZ[i][j] = Vn2->vert[i]->gl[j];
   M1 = long_qbase(M);
   free_mat(M);
   if(M1->rows != fdim)
   {
     printf("error in 'calc_voronoi_stab': form not G-perfect\n");
     exit(3);
   }
   bas = (matrix_TYP **)xmalloc(fdim *sizeof(matrix_TYP *));
   for(i=0;i<fdim;i++)
     bas[i] = vec_to_form(M1->array.SZ[i], Vn2->Gtrred->form, fdim);
   free_mat(M1);

   /*******************************************************************\
   | calculate 'W', the symmetric matrices given by the vertices
   | in Vn1->vert
   \*******************************************************************/
   W = (matrix_TYP **) xmalloc(Vanz *sizeof(matrix_TYP *));
   for(i=0;i<Vanz;i++)
     W[i] = vec_to_form(Vn1->vert[i]->gl, Vn1->Gtrred->form, fdim);

   /*******************************************************************\
   | calculate  an isometry 
   \*******************************************************************/

   /* output for debugging
   put_mat(Vn1->red_inv,NULL,"Vn1->redinv",2);
   put_mat(Vn2->red_inv,NULL,"Vn2->redinv",2);
   put_mat(SV1,NULL,"SV1",2);
   put_mat(SV2,NULL,"SV2",2); */


   X = perfect_normal_isometry(Vn1->red_inv, Vn2->red_inv, SV1, SV2,
          NULL, 0, NULL, W, Vanz, bas, fdim);
   for(i=0;i<fdim;i++)
     free_mat(bas[i]);
   free(bas);
   for(i=0;i<Vanz;i++)
     free_mat(W[i]);
   free(W);
   if(Vn1->SVi == NULL)
     free_mat(SV1);
   free_mat(SV2);

   if(X != NULL)
   {
     M = long_mat_inv(X);
     M1 = mat_mul(Vn1->Ti, M);
     E = mat_mul(M1, Vn2->T);
     free_mat(M1); free_mat(M);
     if(permute == TRUE)
     {
       M = long_mat_inv(E);
       free_mat(E);
       E = M;
       M = NULL;
     }
     free_mat(X);
   }
   else
     E = NULL;
   return(E);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ void calc_voronoi_dir_reps(V, G, Gtr, bifo)
@ voronoi_TYP *V;
@ bravais_TYP *G, *Gtr;
@ matrix_TYP *bifo;
@ 
@  calculates representatives of the orbits of V->stab on the
@  Voronoi-directions given in V->pol->vert.
@  The representatives are stored in V->dir_reps in the following way
@  V->dir_reps->cols is the number of representatives.
@  Each column discribes a representatives
@  The integer in the first row describs the number of the representative,
@  i.e. V->dir_reps->array.SZ[0][i] = k means that V->pol->vert[k]
@  is a representative of the i-th orbit.
@  The entry in the second row is the length of the orbit.
@  The third row is allocated in 'normalizer', where repesentatives
@  V[0],...,V[t] of the G-perfect forms are calculated
@  V->dir_reps->array.SZ[2][i] = k means, that V is a neighbour of
@  V[k] and the Voronoi-direction from V to V[k] is the one
@ given in the first row.
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
calc_voronoi_dir_reps (voronoi_TYP *V, bravais_TYP *G, bravais_TYP *Gtr, matrix_TYP *bifo)
{
   int i,j;
   matrix_TYP **M;
   int fdim;
   int opt[6] = { 1, 0, 0, 0, 0, 0 };
   int orb_no;

   fdim = G->form_no;
   if(V->prime == 0)
    calc_voronoi_basics(V, G, Gtr, 101);
   if(V->stab == NULL)
     calc_voronoi_stab(V, G, Gtr, bifo);
   if(V->pol == NULL)
     calc_voronoi_pol(V, bifo);
   M = (matrix_TYP **) xmalloc(V->pol->vert_no *sizeof(matrix_TYP *));
   for(i=0;i<V->pol->vert_no;i++)
   {
     M[i] = init_mat(1, fdim, "");
     for(j=0;j<fdim;j++)
       M[i]->array.SZ[0][j] = V->pol->vert[i]->v[j];
   }
   V->dir_reps = orbit_representatives(M, V->pol->vert_no, V->linstab, opt, &orb_no, 0);
   for(i=0;i<V->pol->vert_no;i++)
     free_mat(M[i]);
   free(M);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ void calc_voronoi_complete(V, G, Gtr, bifo, prime)
@ voronoi_TYP *V;
@ bravais_TYP *G, *Gtr;
@ matrix_TYP *bifo;
@ int prime;
@
@  applies:
@    calc_voronoi_basics(V, G, Gtr, prime);
@    calc_voronoi_pol(V, bifo);
@    calc_voronoi_good_inv(V, Gtr);
@    calc_voronoi_stab(V, G, Gtr, bifo);
@    calc_voronoi_dir_reps(V, G, Gtr, bifo);
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
calc_voronoi_complete (voronoi_TYP *V, bravais_TYP *G, bravais_TYP *Gtr, matrix_TYP *bifo, int prime)
{
   calc_voronoi_basics(V, G, Gtr, prime);
   calc_voronoi_pol(V, bifo);
   calc_voronoi_good_inv(V, Gtr); 
   calc_voronoi_stab(V, G, Gtr, bifo);
   calc_voronoi_dir_reps(V, G, Gtr, bifo);
}
