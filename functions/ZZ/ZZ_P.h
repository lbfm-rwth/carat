/*
 *  DON'T mangle private and public data, don't include this file if 
 *  _ZZ_H is defined, i.e. ZZ.h has already been included
 */
#if !defined(_ZZ_P_H) && !defined(_ZZ_H)
#define _ZZ_P_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "typedef.h"
#include "matrix.h"
#include "tools.h"
#include "symm.h"
#include "getput.h"

#define  ZZ_PRIVATE

/*
 * this is for the actual implementation. Users won't get much more than
 * the prototype for the function ZZ()
 */

/*{{{  typedef, private types. */

typedef struct {
	int low, hi;
} ZZ_prod_t;

typedef struct ZZ_couple ZZ_couple_t;
typedef struct ZZ_node ZZ_node_t;

struct ZZ_couple {
	ZZ_node_t *he;
	struct {int i, j; } she;
	long factor;
	ZZ_couple_t *elder;
};

struct ZZ_node {
	int number, level, anz_tg;
	int **k_vec;
	ZZ_prod_t *path;
	long index;
	matrix_TYP *U, *U_inv, *el_div;
	ZZ_node_t *next;
	ZZ_couple_t *parent, *child;
	bravais_TYP *group;        /* the groups representation on lattice U */
	bravais_TYP *col_group;    /* transposed of group */
	bravais_TYP *brav;         /* bravais group of col_group */
	bahn **stab_chain;         /* a stabilizer chain for col_group */
	matrix_TYP ***N_orbits;
	int *N_lengths;
	int N_no_orbits;
	voronoi_TYP **perfect;
	int perfect_no;
};

/*--------------------------------------------------------------------*\
| r    = #generators                  k = #primes                      |
| p[i] = i-th prime number         s[i] = #constituents for i-th prime |
| n[i][j]        = #rows of j-th constituent for i-th prime            |
| Delta[i][j][k] = j-th constituent for i-th prime k-th generator      |
| Endo [i][j] =   endomorphisms for Delta[i][j] and all generators     |
| EnCo [i][j] = # endomorphisms for Delta[i][j]  "   "      "          |
| VK[i] = Vielfachheitenvektor der Konstituenten zur i-ten Primzahl    |
\*--------------------------------------------------------------------*/
typedef struct {
	int k;
	int *s, *p;
	matrix_TYP ****Delta;
} ZZ_prime_constituents_t;

typedef struct {
	ZZ_prime_constituents_t p_consts;
	int N, r;
	int **n;
	ZZ_prod_t **EnCo;
	matrix_TYP **DELTA_M, **DELTA;
	matrix_TYP **epi_base, *epi;
	matrix_TYP ****Endo;
	int **VK;
} ZZ_data_t;

typedef struct {
	ZZ_node_t *root, *last;
	int node_count;
} ZZ_tree_t;

/*{{{  global variables */
extern boolean QUIET;
extern boolean TEMPORAER;
extern boolean SHORTLIST;
extern boolean NURUMF;
extern boolean U_option;
extern boolean G_option;
extern boolean LLLREDUCED;
extern int COUNTER;
extern int NUMBER;
extern int ABBRUCH;
extern int ZCLASS;
extern int SUBDIRECT;
extern int LEVEL;
extern FILE *ZZ_temp;
extern FILE *ZZ_list;
extern int MAT_ALLOC;
extern int constituents;
extern int verbose;

extern int IDEM_NO;
/*}}}  */

#if defined(__STDC__)
#define _ZZ_P_PROTO_( args ) args
#else
#define _ZZ_P_PROTO_( args ) ()
#endif

extern void ZZ_transpose_array _ZZ_P_PROTO_((int **array, int size));

extern void ZZ_intern _ZZ_P_PROTO_((matrix_TYP * Gram,
				    ZZ_data_t * data,
				    ZZ_tree_t * tree)) ;
extern void ZZ _ZZ_P_PROTO_((bravais_TYP * group, 
			     matrix_TYP * gram, 
			     int *divisors, 
			     char *options,
			     FILE *putputfile)) ;

#endif /* _ZZ_P_H */

