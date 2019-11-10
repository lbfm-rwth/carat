#include "typedef.h"
#include "tools.h"
#include "matrix.h"
#include "voronoi.h"
#include "ZZ_P.h"
#include "ZZ_gen_vs_P.h"

/*
 *  Determine the image of the generators in Mn(GF(p),n). 
 */

matrix_TYP **ZZ_mod_gen (matrix_TYP **gen, int num)
{
	int i, j, k, dim;
	matrix_TYP **p_gen;

	dim = gen[0]->rows;
	p_gen = (matrix_TYP **) malloc (num * sizeof (matrix_TYP *));
	for (k = 0; k < num; k++) {
		p_gen[k] = init_mat (dim, dim, "p");
		for (i = 0; i < dim; i++) {
			for (j = 0; j < dim; j++) {
				p_gen[k]->array.SZ[i][j] = 
					gen[k]->array.SZ[i][j] % act_prime;
				if (p_gen[k]->array.SZ[i][j] < 0) {
					p_gen[k]->array.SZ[i][j] += act_prime;
				}
			}
		}
	}
	return (p_gen);
}


/*
 * Determine the rank and a basis for the image of a vector 
 * under the operation of an order (given by its Z-basis).\
 */

matrix_TYP *ZZ_vec_bahn (int *vec, matrix_TYP **gen, int num)
{
	matrix_TYP *basis;
	int **B, **T, *v, new_flag;
	int i, j, k, dim, act, pos, rank, flag;
	
	dim = gen[0]->rows;
	
	basis = init_mat ((num + 1) * dim, dim, "p");
	basis->prime = gen[0]->prime;
	B = basis->array.SZ;
	memcpy (B[0], vec, dim * sizeof (int));
	rank = 1;
	pos = 1;
	flag = 1;
	while (flag) {
		for (act = 0; act < rank; act++) {
			v = B[act];
			for (i = 0; i < num; i++) {
				if (save_null_mat (gen[i]))
					continue;
				T = gen[i]->array.SZ;
				new_flag = 0;
				for (j = 0; j < dim; j++) {
					B[pos][j] = P (v[0], T[0][j]);
					for (k = 1; k < dim; k++) {
						B[pos][j] = S(B[pos][j], P(v[k], T[k][j]));
					}
					if (B[pos][j])
						new_flag++;
				}
				if (new_flag)
					pos++;
			}
		}
		basis->rows = pos;
		i = p_gauss (basis);
		flag = ((i == rank) || (i == dim)) ? 0 : 1;
		rank = i;
		pos = i;
	}
	for (i = rank; i < (num + 1) * dim; i++) {
		free (basis->array.SZ[i]);
	}
	basis->array.SZ = (int **) realloc (basis->array.SZ, rank * sizeof (int *));
	basis->rows = rank;
	
	Check_mat (basis);
	return (basis);
}


/*
 * Generate the vector space of given dimension over the
 * prim field GF(p) up to a skalar factor (the leading
 * coordinate unequal to zero is one).
 */

matrix_TYP *ZZ_gen_vs (int prim, int dim)
{
	int act, i, j, n;
	matrix_TYP *VS;
	int **M;
	

	/* calculate number of vectors */
	n = 1;
	j = prim;
	for (i = 1; i < dim; i++) {
		n += j;
		j *= prim;
	}
	
	/* alloc space and generated vector space */
	
	VS = init_mat (n, dim, "");
	M = VS->array.SZ;
	i = 1;
	act = dim - 1;
	M[0][act] = 1;
	while (i < n) {
		j = dim - 1;
		while ((j > 0) && ((M[i][j] = (M[i - 1][j] + 1) % prim) == 0))
			j--;
		if (j <= act) {
			M[i][act--] = 0;
			M[i][act] = 1;
		} else {
			while (--j >= act)
				M[i][j] = M[i - 1][j];
		}
		i++;
	}
	return (VS);
}
