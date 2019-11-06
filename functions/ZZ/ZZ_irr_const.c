#include "typedef.h"
#include "matrix.h"
#include "voronoi.h"
#include "tools.h"
#include "ZZ_P.h"
#include "ZZ_gen_vs_P.h"

#define TOPLEVEL -1

#define TRACE(level, foo)			\
if (level <= TOPLEVEL) {			\
	foo;					\
}

/*  Calculate the irreducible constituents of an order as input for the
 *  centering algorithm (to determine the invariant lattices) and for
 *  determine the radical.  The order is supposed to be normalized,
 *  i.e. the standard maximal order Mn(R,dim) is an over_order of the
 *  given one.  
 *  
 *  first a help-programm: Do the calculation for a matrix- set of
 *  generators and use this recursivly.  
 */

/*
 *                                                                    
 *  input:                                                              
 *  matrix_TYP **generators:                                          
 *     die Generatoren, fuer jeden Generator eine Matrix              
 *  int num_gens: Anzahl der Generatoren                              
 *  int *num_irr_const: die Anzahl der irrduziblen Konstituenten.              
 *                                                                    
 *   output:                                                          
 *   matrix_TYP **:    die irreduziblen Konstituenten mod p fuer      
 *                     jeden Generator. Zuerst die Matrizen aller     
 *                     Generatoren fuer den ersten Konstituenten,     
 *                     dann fuer den zweiten, usf.                    
 *                                                                    
 */

static matrix_TYP **irr_help(matrix_TYP ** generators,
			      int num_gens,
			      int *num_irr_const,
			      int rec_cnt)
{
	matrix_TYP *help, *mat, **upset1, **upset2, **gen1, **gen2;
	matrix_TYP **r_gen, **p_gen, *VS;
	int dim, rg, notfound, i, j, act, cnt;
	int **V, **T1, **T2;

	/* p_gen enthaelt die Generatoren mod p */
	p_gen = ZZ_mod_gen(generators, num_gens);
	for (i = 0; i < num_gens; i++) {
		TRACE(1, fprintf(stderr, "generator nr %d\n", i));
		TRACE(1, fput_mat(stderr, p_gen[i], "generator mod p", 0));
	}
	if (p_gen[0]->rows == 1) { /* really nothing to do */
		*num_irr_const = 1;
		return p_gen;
	}
	dim = p_gen[0]->rows;
	VS = ZZ_gen_vs(act_prime, dim);
	TRACE(2, fput_mat(stderr, VS, "Vektorraum", 0));
	for (notfound = 1, i = 0; (i < VS->rows) && notfound; i++) {
		/* help enthaelt eine Basis der VR, der von der Bahn
		 * aufgespannt wird 
		 */
		help = ZZ_vec_bahn(VS->array.SZ[i], p_gen, num_gens);
		TRACE(2, fput_mat(stderr, help, "Bahn", 0));
		notfound = help->rows == dim;
		if (notfound) {
			free_mat(help);
		}
	}
	free_mat(VS);
	if (notfound) {
		*num_irr_const = 1;
		TRACE(0,fprintf(stderr, "Kein Basiswechsel\n"));
		return p_gen;
	}
	r_gen=(matrix_TYP **)malloc(num_gens * sizeof(matrix_TYP *));
	TRACE(1,fprintf(stderr,
		      "Found submodule of dimension %d, recursion count %d.\n",
			help->rows, rec_cnt));
	/*  Submodule found: begin the recursion.
	 *  Start with filling up help to a basis
	 */
	mat = init_mat(dim, dim, "p");
	mat->prime = act_prime;
	rg = help->rows;
	act = 0;
	TRACE(1, fput_mat(stderr, help, "Invarianter Teilraum", 0));
	for (i = 0; i < help->rows; i++) {
		memcpy(mat->array.SZ[i],
			help->array.SZ[i], dim * sizeof(int));
		while(mat->array.SZ[i][act] == 0) {
			mat->array.SZ[rg][act] = 1;
			rg++;
			act++;
		}
		act++;
	}
	for (; rg < mat->rows; rg++, act++) {
		mat->array.SZ[rg][act] = 1;
	}
	TRACE(2, fput_mat(stderr, mat, "Basis?", 0));
	rg = help->rows;
	free_mat(help);
	help = mat_inv(mat);
	TRACE(0,fprintf(stderr, "Recursion: %d", rec_cnt));
	TRACE(0,fput_mat(stderr, mat, "Basiswechsel", 0));

	/*  Do the projection to the submoduls
	 */
	upset1 = (matrix_TYP **) malloc(num_gens * sizeof(matrix_TYP *));
	upset2 = (matrix_TYP **) malloc(num_gens * sizeof(matrix_TYP *));
	for (i = 0; i < num_gens; i++) {
		/* Basiswechsel fuer die Generatoren durchfuehren
		 * diese operieren dann auf den Standardbasisvektoren
		 * so dass gen_vc wieder benutzt werden kann  
		 * diese sollten jetzt eine Blockmatrix darstellen.
		 */
		r_gen[i] = mat_kon(mat, p_gen[i], help);
		TRACE(2,fprintf(stderr,
				"Transformierter Generator %d mod p, rg: %d\n",
				i, rg));
		TRACE(2,fput_mat(stderr, r_gen[i], "Generator transformiert",
				  0));
		upset1[i] = init_mat(rg, rg, "ik");
		upset2[i] = init_mat(dim - rg, dim - rg, "ik");
		upset1[i]->prime = act_prime;
		upset2[i]->prime = act_prime;
		V = r_gen[i]->array.SZ;
		T1 = upset1[i]->array.SZ;
		T2 = upset2[i]->array.SZ;
		for (j = 0; j < dim; j++) {
			if (j < rg) {
				/* Generatoren fuer den
				   invarianten Teilraum */
				memcpy(T1[j], V[j],rg * sizeof(int));
			} else {
				/* Generatoren fuer den Rest */
				memcpy(T2[j - rg], &V[j][rg],
					(dim - rg) * sizeof(int));
			}
		}
		/* die neuen Generatoren sind in upset_i */
		free_mat(r_gen[i]);
		free_mat(p_gen[i]);
	}
	free_mat(mat);
	free_mat(help);
	free(r_gen);
	free(p_gen);
	/*----------------------------------------------------*\
	  | Create the recursion-input                           |
	  \*----------------------------------------------------*/
	gen1 = upset1;
	gen2 = upset2;
	TRACE(0, fprintf(stderr, "Teilraum\n"));
	upset1 = irr_help(upset1, num_gens, &i, rec_cnt + 1);
	TRACE(0, fprintf(stderr, "Faktorraum\n"));
	upset2 = irr_help(upset2, num_gens, &j, rec_cnt + 1);
	for (cnt = 0; cnt < num_gens; cnt++) {
		free_mat(gen1[cnt]);
		free_mat(gen2[cnt]);
	}
	free(gen1);
	free(gen2);
	*num_irr_const = i + j;
	p_gen = (matrix_TYP **) malloc(*num_irr_const * num_gens * sizeof(matrix_TYP *));
	memcpy(p_gen, upset1, i * num_gens * sizeof(matrix_TYP *));
	memcpy(p_gen + num_gens * i, upset2,
	       j * num_gens * sizeof(matrix_TYP *));
	free(upset1);
	free(upset2);
	return p_gen;
}

/***********************************************************************/
/*                                                                     */
/* input:                                                              */
/*   matrix_TYP **generators:                                       */
/*      die Generatoren, fuer jeden Generator eine Matrix              */
/*   int num_gens: Anzahl der Generatoren                              */
/*   int *num_irr_const: die Anzahl der irrduziblen Konstituenten.     */
/*   int p: die Primzahl, modulo derer die Konstituenten berechnet     */
/*          werden                                                     */
/*                                                                     */
/*    output:                                                          */
/*    matrix_TYP **: die irreduziblen Konstituenten mod p fuer      */
/*                      jeden Generator. Zuerst die Matrizen aller     */
/*                      Generatoren fuer den ersten Konstituenten,     */
/*                      dann fuer den zweiten, usf.                    */
/*                                                                     */
/***********************************************************************/

matrix_TYP **ZZ_irr_const(matrix_TYP **generators,
			   int num_gens, int p, int *num_irr_const)
{
	matrix_TYP **irrconst;
	int i, j;

	init_prime(p);
	irrconst = irr_help(generators, num_gens, num_irr_const, 0);
#ifdef DEBUG
	fprintf(stderr, "Primzahl: %d\n",p );
#endif
	for (i = 0; i < *num_irr_const; i++) {
		for (j = 0; j < num_gens; j++) {
			TRACE(1, fprintf(stderr, "Konstituent: %d\n", i));
			TRACE(1, fprintf(stderr, "Generator: %d\n", j));
			TRACE(1, fput_mat(stderr, irrconst[i*num_gens + j], "                          :", 2));
		}
	}
	return irrconst;
}
