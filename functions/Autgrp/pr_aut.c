/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  pr_aut.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

#include "typedef.h"
#include "matrix.h"
#include "reduction.h"
#include "symm.h"
#include "autgrp.h"
#include "types.h"
#include "datei.h"

static int normal_option;
static int perp_no;
static int ***perp;
static int perpdim;
static int ***perpbase;
static int ***perpprod;
static int *perpvec;


/*************************************************************************\
| The functions 'pr_aut' calculates generators and the order of the group
| G := {g in GL_n(Z) | g * Fo[i] * g^{tr} = Fo[i], 1<= i<= Foanz}
| returned via a pointer to 'bravais_TYP'.
| This function applues a pair reduction to Fo[0] before calculating
| the group
|
| The arguments of pr_aut are:
| matrix_TYP	**Fo:		a set of n times n matrices,
|				the first must be positiv definite
| int	 	Foanz:		the number of the matrices given in 'Fo'.
| matrix_TYP	**Erz:		if already element of G are known,
|				they can be used for calculating generators
|				for the whole group.
|				The matrices of known elements can be given 
|				to the function by the pointer 'Erz'.
| int		Erzanz:		The number of matrices given in 'Erz"
| int		*options:	see below.
|
| options is a pointer to integer (of length 6)
| The possible options are encoded in the following way:
| options[0]:	The depth, up to wich scalar product combinations
|		shall be calculated. The value should be small.
|		options[0] > 0 should be used only, if the automorphismn
|		group is expected to be small (with respect to the number
|		of shortest vectors).
| options[1]:	The n-point stabiliser with respect to different basis
|		will be calculated.
| options[2]:	If options[2] = 1, additional output is written to the
|               file AUTO.tmp
| options[3]:   If options[3] = 1, Bacher polynomials are used. 
|		If options[3] = 2, Bacher polynomial are used up to a deepth
|		                   specified in options[4].
|		If options[3] = 3, Bacher polynomials are used, using 
|		                   combinations of vectors having the scalar
|		                   product specified in options[5]
|		options[3] = 4 is the combination of options[3] = 2 and
|                              options[3] = 3.
| options[4]:	A natural number number  or zero (if options[3] = 2 or 4)
| options[5]:	An integral number (if options[3] = 3 or 4)
|
|	It is possible to use NULL for options,
|	in this case option is assumed to be [0,0,0,0,0,0]
\*************************************************************************/


bravais_TYP *pr_aut(Fo, Foanz, Erz, Erzanz, options)
matrix_TYP **Fo, **Erz;
int Foanz, Erzanz, *options;
{
   matrix_TYP **F, **E, *SV, *T, *Titr, *Ttr, *w;
   bravais_TYP *G, *G1;
   int n, anz, i;

   if((F = (matrix_TYP **)malloc(Foanz *sizeof(matrix_TYP *))) == NULL) 
   {
     printf("malloc of F in 'pr_aut' failed\n");
     exit(2);
   }
   E = NULL;
   if(Erzanz != 0)
     if((E = (matrix_TYP **)malloc(Foanz *sizeof(matrix_TYP *))) == NULL) 
     {
       printf("malloc of E in 'pr_aut' failed\n");
       exit(2);
     } 
   T = init_mat(Fo[0]->rows, Fo[0]->cols, "");
   n = Fo[0]->cols;
   for(i=0;i<n;i++)
     T->array.SZ[i][i] = 1;
   F[0] = pair_red(Fo[0], T);
   Ttr = tr_pose(T);
   Titr = mat_inv(Ttr);
   for(i=1;i<Foanz;i++)
       F[i] = scal_pr(T, Fo[i], 1);
   SV = short_vectors(F[0], F[0]->array.SZ[n-1][n-1], 0, 0, 0, &anz);
   for(i=0;i<Erzanz;i++)
   {
     w = mat_mul(Titr, Erz[i]);
     E[i] = mat_mul(w, Ttr);
     free_mat(w);
   }
   G1 = autgrp(F, Foanz, SV, E, Erzanz, options);
   if ((G = (bravais_TYP *)malloc(sizeof(bravais_TYP))) == 0)
   {
     printf("malloc of 'G' in 'pr_aut' failed\n");
     exit(2);
   }
   G->gen_no = G1->gen_no;
   G->form_no = G->normal_no = G->cen_no = G->zentr_no = 0;
   for(i=0;i<100;i++)
     G->divisors[i] = G1->divisors[i];
   G->order = G1->order;
   G->dim = G1->gen[0]->cols;
   if((G->gen = (matrix_TYP **)malloc(G->gen_no *sizeof(matrix_TYP *))) == NULL)
   {
     printf("mallocof 'G->gen' in 'pr_aut' failed\n");
     exit(2);
   }
   for(i=0;i<G->gen_no;i++)
   {
     w = mat_mul(Ttr, G1->gen[i]);
     G->gen[i] = mat_mul(w, Titr);
     free_mat(w);
   }
   free_mat(T);
   free_mat(Ttr);
   free_mat(Titr); free_mat(SV);
   for(i=0;i<Erzanz;i++)
     free_mat(E[i]);
   if(Erzanz != 0)
     free(E);
   for(i=0;i<Foanz;i++)
     free_mat(F[i]);
   free(F);
   free_bravais(G1);
   return(G);
}



/*********************************************
bravais_TYP *perfect_normal_autgrp(Fo, SV, Erz, Erzanz, options, P, Panz, Pbase, Pdim)
matrix_TYP *Fo, **Erz, *SV, **P, **Pbase;
int Erzanz, *options, Panz, Pdim;
*********************************************/
