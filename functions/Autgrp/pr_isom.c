/*****	Main program for the isometry program ISOM	*****/
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: pr_isom.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

#include "typedef.h"
#include "types.h"
#include "matrix.h"

static int normal_option;
static int perp_no;
static int ***perp;
static int perpdim;
static int ***perpbase;
static int ***perpprod;
static int *perpvec;
/*
@-------------------------------------------------------------------------
@ matrix_TYP *pr_isom(F1, F2, Fanz, Erz, Erzanz, options)
@ matrix_TYP **F1, **F2, **Erz;
@ int Fanz, Erzanz, *options;
@
@ The function 'isometry' calculates a matrix X such that
@    X * F1[i] *X^{tr} = F2[i] for 1<= i<= Foanz
@ returned via a pointer to 'matrix_TYP'.
@ If no such matrix exists the functions returns NULL
@ 'pr_aut' applies a pair_reduction to F1[0] and then
@ uses then the function isometry to calculate an isometry
@ for the reduced form
@
@ The arguments of isometry are:
@ matrix_TYP	**F1:		a set of n times n matrices,
@				the first must be positiv definite
@ matrix_TYP	**F2:		a set of n times n matrices,
@				the first must be positiv definite
@ int	 	Fanz:		the number of the matrices given in 'F1'.
@ int		*options:	see below.
@ matrix_TYP	**Erz:		if already elements with g^{tr}F1[i]g = F2[i]
@                               are known, the can be used for calculating
@                               the isometry.
@				The matrices of known elements can be given 
@				to the function by the pointer 'Erz'.
@ int		Erzanz:		The number of matrices given in 'Erz"
@ int		*options:	see below.
@
@ options is a pointer to integer (of length 6)
@ The possible options are encoded in the following way:
@ options[0]:	The depth, up to wich scalar product combinations
@		shall be calculated. The value should be small.
@		options[0] > 0 should be used only, if the automorphismn
@		group is expected to be small (with respect to the number
@		of shortest vectors).
@ options[1]:	The n-point stabiliser with respect to different basis
@		will be calculated.
@ options[2]:	If options[2] = 1, additional output is written to the
@               file AUTO.tmp
@ options[3]:   If options[3] = 1, Bacher polynomials are used. 
@		If options[3] = 2, Bacher polynomial are used up to a deepth
@		                   specified in options[4].
@		If options[3] = 3, Bacher polynomials are used, using 
@		                   combinations of vectors having the scalar
@		                   product specified in options[5]
@		options[3] = 4 is the combination of options[3] = 2 and
@                              options[3] = 3.
@ options[4]:	A natural number number  or zero (if options[3] = 2 or 4)
@ options[5]:	An integral number (if options[3] = 3 or 4)
@
@	It is possible to use NULL for options,
@	in this case option is assumed to be [0,0,0,0,0,0]
@-------------------------------------------------------------------------
*/

matrix_TYP *pr_isom(F1, F2, Fanz, Erz, Erzanz, options)
matrix_TYP **F1, **F2, **Erz;
int Fanz, Erzanz, *options;
{
  matrix_TYP *T, *X, *X1, **F, *SV1, *SV2;
  int i,j,k;
  int n, anz, max;

  extern matrix_TYP *isometry();
  extern matrix_TYP *short_vectors();
  extern matrix_TYP *mat_mul();
  extern matrix_TYP *scal_pr();
  extern matrix_TYP *pair_red();
  extern matrix_TYP *init_mat();

  if((F = (matrix_TYP **) malloc(Fanz *sizeof(matrix_TYP *))) == NULL)
  {
    printf("malloc of 'F' in 'pr_isom' failed\n");
    exit(2);
  }
  n = F1[0]->cols;
  T = init_mat(n,n,"");
  for(i=0;i<n;i++)
   T->array.SZ[i][i] = 1;
  F[0] = pair_red(F1[0], T);
  max = F[0]->array.SZ[n-1][n-1];
  for(i=1;i<Fanz;i++)
    F[i] = scal_pr(T, F1[i], 1);
  SV1 = short_vectors(F[0], max, 0, 0, 0, &anz);
  SV2 = short_vectors(F2[0], max, 0, 0, 0, &anz);
  X1 = isometry(F, F2, Fanz, SV1, SV2, Erz, Erzanz, options);
  for(i=0;i<Fanz;i++)
     free_mat(F[i]);
  free(F);
  free_mat(SV1); free_mat(SV2);
  if(X1 == NULL)
  {
    free_mat(T); 
    return(NULL);
  }
  X = mat_mul(X1, T);
  free_mat(X1); free_mat(T);
  return(X);
}


/*********************************************************
matrix_TYP *perfect_normal_isometry(F1, F2, SV1, SV2, Erz, Erzanz, options, P, Panz, Pbase, Pdim)
matrix_TYP *F1, *F2, *SV1, *SV2, **Erz, **P, **Pbase;
int Erzanz, *options, Panz, Pdim;
********************************************************/
