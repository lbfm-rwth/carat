/*****	This file includes some routines for 
	the preprocessing of the algorithm	*****/
#include "typedef.h"
#include "utils.h"
#include "types.h"

/************************************************************************\
|	removes those vectors from V.v, which 
|	have other norm combinations than the base-vectors
\************************************************************************/
void 
checkvecs (veclist *V, invar F, veclist norm)
{
	int	i, j, k, dim, *normvec, *Vvi;

	dim = F.dim;
	normvec = (int*)xmalloc(norm.dim * sizeof(int));
	j = 0;
	for (i = 1; i <= V->n; ++i)
	{
		Vvi = V->v[i];
		for (k = 0; k < F.n; ++k)
		{
			Vvi[dim+k] = scp(Vvi, F.A[k], Vvi, dim);
			normvec[k] = Vvi[dim+k];
		}
		if (abs(numberof(normvec, norm)) > dim)
/* the vector V->V[i] has a different norm combination from those in the list
   norm.v and is therefore deleted */
		{
			free(Vvi);
			++j;
		}
		else
			V->v[i-j] = Vvi;
	}
	V->n -= j;
	free(normvec);
}

/************************************************************************\
|	checks, whether g fixes the invariant forms
\************************************************************************/
int 
checkgen (int **g, invar F)
{
	int	i, j, k, fix, **FAk;

	fix = 1;
	for (k = 0; k < F.n  &&  fix == 1; ++k)
	{
		FAk = F.A[k];
		for (i = 0; i < F.dim  &&  fix == 1; ++i)
		{
			for (j = 0; j <= i  &&  fix == 1; ++j)
			{
				if (scp(g[i], FAk, g[j], F.dim) != FAk[i][j])
					fix = 0;
			}
		}
	}
	return(fix);
}

/************************************************************************\
| a permutation per := fp.per for the order of the basis-vectors is chosen 
| such that in every step the number of possible continuations is minimal, 
| for j from per[i] to per[dim-1] the value f[i][j] in the fingerprint f is 
| the number of vectors, which have the same scalar product with the  
| basis-vectors per[0]...per[i-1] as the basis-vector j and the same length as 
| this vector with respect to all  invariant forms
\************************************************************************/
void 
fingerprint (fpstruct *fp, invar F, veclist V)
{
	int	i, j, k, dim, min, tmp, **f, *Vvj;

	dim = F.dim;
	f = (int**)xmalloc(dim * sizeof(int*));
	for (i = 0; i < dim; ++i)
	{
		f[i] = (int*)xmalloc(dim * sizeof(int));
		for (j = 0; j < dim; ++j)
			f[i][j] = 0;
	}
	for (i = 0; i < dim; ++i)
		fp->per[i] = i;
/* the first row of the fingerprint has as entry nr. i the number of
   vectors, which have the same length as the i-th basis-vector with
   respect to every invariant form */
	for (j = 1; j <= V.n; ++j)
	{
		Vvj = V.v[j];
		for (i = 0; i < dim; ++i)
		{
			for (k = 0; k < F.n  &&  Vvj[dim+k] == F.A[k][i][i]; ++k);
			if (k == F.n)
				f[0][i] += 2;
		}
	}
	for (i = 0; i < dim-1; ++i)
	{
/* a minimal entry != 0 in the i-th row is chosen */
		min = i;
		for (j = i+1; j < dim; ++j)
		{
			if (f[i][fp->per[j]] < f[i][fp->per[min]])
				min = j;
		}
		tmp = fp->per[i];
		fp->per[i] = fp->per[min];
		fp->per[min] = tmp;
/* the column below the minimal entry is set to 0 */
		for (j = i+1; j < dim; ++j)
			f[j][fp->per[i]] = 0;
/* compute the row i+1 of the fingerprint */
		for (j = i+1; j < dim; ++j)
			f[i+1][fp->per[j]] = possible(F, V, fp->per, i, fp->per[j]);
	}
/* only the diagonal of f will be needed later */
	for (i = 0; i < dim; ++i)
	{
		fp->diag[i] = f[i][fp->per[i]];
		free(f[i]);
	}
	free(f);
}

/************************************************************************\
|	returns the number of vectors,
|	which have the same scalar products with the 
|	basis-vectors nr. per[0]...per[I] as the 
|	basis-vector nr. J and the same length 
|	as this vector with respect to all invariant forms	
\************************************************************************/
int 
possible (invar F, veclist V, int *per, int I, int J)
{
	int	i, j, k, dim, count, *Vvj;

	dim = F.dim;
	count = 0;
	for (j = 1; j <= V.n; ++j)
	{
		Vvj = V.v[j];
		i = I+1;
/* check the length of the vector */
		for (k = 0; k < F.n  &&  i > I  &&  Vvj[dim+k] == F.A[k][J][J]; ++k)
/* check the scalar products with the basis-vectors */
			for (i = 0; i <= I  &&  F.v[k][j][per[i]] == F.A[k][per[i]][J]; ++i);
		if (k == F.n  &&  i > I)
			++count;
/* the same for the negative vector */
		i = I+1;
		for (k = 0; k < F.n  &&  i > I  &&  Vvj[dim+k] == F.A[k][J][J]; ++k)
			for (i = 0; i <= I  &&  F.v[k][j][per[i]] == -F.A[k][per[i]][J]; ++i);
		if (k == F.n  &&  i > I)
			++count;
	}
	return(count);
}

/************************************************************************\
|	calculates the scalar products of the vector w with the base
|	vectors v[b[I]] down to v[b[I-dep+1]] with respect to
|	all invariant forms and puts them on scpvec
\************************************************************************/
void 
scpvector (int *scpvec, int *w, int *b, int I, int dep, invar F)
{
	int	i, j, dim, bi;

	dim = F.dim;
	for (i = I; i >= 0  &&  i > I-dep; --i)
	{
		if ((bi = b[i]) > 0)
		{
			for (j = 0; j < F.n; ++j)
				scpvec[j*dep + I-i] = sscp(w, F.v[j][bi], dim);
		}
		else
		{
			for (j = 0; j < F.n; ++j)
				scpvec[j*dep + I-i] = -sscp(w, F.v[j][-bi], dim);
		}
	}
}

/************************************************************************\
|	computes the list of scalar product 
|	combinations of the vectors in V.v with
|	the basis-vectors in b
\************************************************************************/
void 
scpvecs (veclist *list, int ***vec, int I, int *b, int dep, veclist V, invar F)
{
	int	i, j, dim, len, *scpvec, nr, sign, *tmp;
	int	*Vvj, *vecnr, *listvn, *vecn, **listv;

	dim = F.dim;
	len = F.n * dep;
/* scpvec is the vector for the scalar product combination */
	scpvec = (int*)xmalloc(len * sizeof(int));
/* the first vector in the list is the 0-vector and is not counted */
	list->v = (int**)xmalloc(1 * sizeof(int*));
	list->dim = len;
	list->len = len;
	list->v[0] = (int*)xmalloc(len * sizeof(int));
	for (i = 0; i < len; ++i)
		list->v[0][i] = 0;
	list->n = 0;
	*vec = (int**)xmalloc(1 * sizeof(int*));
	(*vec)[0] = (int*)xmalloc(dim * sizeof(int));
	for (i = 0; i < dim; ++i)
		(*vec)[0][i] = 0;
	for (j = 1; j <= V.n; ++j)
	{
		Vvj = V.v[j];
		for (i = 0; i < len; ++i)
			scpvec[i] = 0;
		scpvector(scpvec, Vvj, b, I, dep, F);
		for (i = 0; i < len  &&  scpvec[i] == 0; ++i);
		if (i == len)
			continue;

        nr = numberof(scpvec, *list);
        sign = nr > 0 ? 1 : -1;
        nr = abs(nr);
/* scpvec is already in list */
		if (nr <= list->n  &&  nr > 0)
		{
			vecnr = (*vec)[nr];
			for (i = 0; i < dim; ++i)
				vecnr[i] += sign * Vvj[i];
		}
/* scpvec is a new scalar product combination */
		else if (nr >= list->n+1)
		{
			++list->n;
			list->v = (int**)xrealloc(list->v, (list->n+1) * sizeof(int*));
			list->v[list->n] = (int*)xmalloc(len * sizeof(int));
/* numberof changes the sign of scpvec if the first nonzero entry is < 0,
   hence this is correct */
			listvn = list->v[list->n];
			for (i = 0; i < len; ++i)
				listvn[i] = scpvec[i];
			*vec = (int**)xrealloc(*vec, (list->n+1) * sizeof(int*));
			(*vec)[list->n] = (int*)xmalloc(dim * sizeof(int));
			vecn = (*vec)[list->n];
			for (i = 0; i < dim; ++i)
				vecn[i] = sign * Vvj[i];
/* shuffle the new vector to the correct position, this should be quick enough,
   since the length of the list should not exceed some hundreds */
			listv = list->v;
			for (i = list->n; i > nr+1-list->n; --i)
			{
				tmp = listv[i];
				listv[i] = listv[i-1];
				listv[i-1] = tmp;
				tmp = (*vec)[i];
				(*vec)[i] = (*vec)[i-1];
				(*vec)[i-1] = tmp;
			}
		}
	}
	free(scpvec);
}

/************************************************************************\
|	computes a basis b for the lattice generated by the vectors in v,
|	puts a transformation matrix on T, i.e. b = T*v,
|	uses LLL-reduction with respect to F
\************************************************************************/
void 
base (scpcomb *com, int ***b, int **v, int **F, int dim)
{
	int	i, j, k, nv, rank, max, **Fv, **f, **tr, *perm, *norm, tmp;
	int	*vi, *Fvi, *comtr, *fj;

	nv = com->list.n + 1;
	max = 1;
/* get the maximal entry in the vector sums */
	for (i = 0; i < nv; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			if (abs(v[i][j]) > max)
				max = abs(v[i][j]);
		}
	}
/* Fv is the list of products of F with the vectors in v, scalar products
   with respect to F can be calculated as standard scalar products of
   vectors in v with vectors in Fv */
	Fv = (int**)xmalloc(nv * sizeof(int*));
/* the product of the maximal entry in the vector sums with the maximal entry
   in the product of the Gram-matrices with the vector sums should not exceed
   MAXNORM to avoid overflow */
	max = MAXNORM / max;
	for (i = 0; i < nv; ++i)
	{
		Fv[i] = (int*)xmalloc(dim * sizeof(int));
		Fvi = Fv[i];
		for (j = 0; j < dim; ++j)
		{
			Fvi[j] = sscp(F[j], v[i], dim);
			if (abs(Fvi[j]) > max)
/* an entry in F.v[i] is too large */
			{
				fprintf(stderr, "Error: Found entry %d in vector sum.\nTo avoid overflow, the entries should not exceed %d.\nTry without option -D or use -D with higher value.\n", Fvi[j], max);
				exit (3);
			}
		}
	}
/* b is the basis for the lattice */
	*b = (int**)xmalloc(1 * sizeof(int*));
	(*b)[0] = (int*)xmalloc(dim * sizeof(int));
/* com->trans is the transformation matrix */
	com->trans = (int**)xmalloc(1 * sizeof(int*));
	com->trans[0] = (int*)xmalloc(nv * sizeof(int));
/* f is the Gram-matrix for the basis b with respect to the form F */
	f = (int**)xmalloc(1 * sizeof(int*));
	f[0] = (int*)xmalloc(1 * sizeof(int));
/* tr is the transformation matrix for the LLL-reduced lattice */
	tr = (int**)xmalloc(1 * sizeof(int*));
	tr[0] = (int*)xmalloc(1 * sizeof(int));
/* perm is the new order in which the vectors in v are added to the lattice
   generated by the preceding vectors */
	perm = (int*)xmalloc(nv * sizeof(int));
	norm = (int*)xmalloc(nv * sizeof(int));
/* bubble sort with respect to the norm of the vectors in v with respect to the
   Gram-matrix F, which is assumed to be positive definite, 
   this should stabilize the LLL-reduction, 
   the bubble sort should be sufficient since the number of vectors here is 
   assumed to be rather small (at most some hundreds) */
	for (i = 0; i < nv; ++i)
	{
		norm[i] = sscp(v[i], Fv[i], dim);
		if (norm[i] > MAXNORM)
		{
			exit (3);
		}
		perm[i] = i;
	}
	i = 0;
	while (i < nv-1)
	{
		if (norm[perm[i]] > norm[perm[i+1]])
		{
			tmp = perm[i];
			perm[i] = perm[i+1];
			perm[i+1] = tmp;
			if (i > 0)
				--i;
		}
		else
			++i;
	}
	free(norm);
/* jump over possible 0-vectors */
	i = 0;
	for (j = 0; j < dim  &&  v[perm[i]][j] == 0; ++j);
	while (j == dim  &&  i < nv)
	{
		++i;
		if (i < nv)
			for (j = 0; j < dim  &&  v[perm[i]][j] == 0; ++j);
	}
	rank = 0;
	while (i < nv)
	{
/* v[perm[i]] is the candidate to enlarge the lattice generated by the 
   base vectors in b */
		vi = v[perm[i]];
		Fvi = Fv[perm[i]];
		for (j = 0; j < rank; ++j)
		{
			f[j][rank] = sscp((*b)[j], Fvi, dim);
			f[rank][j] = f[j][rank];
		}
		f[rank][rank] = sscp(vi, Fvi, dim);
		if (lll(f, tr, rank+1) > rank)
/* the rank of the lattice generated by v[perm[i]] and the vectors in b is 
   rank+1, i.e. one higher than without v[perm[i]] */
		{
			++rank;
			for (j = 0; j < dim; ++j)
				(*b)[rank-1][j] = vi[j];
			comtr = com->trans[rank-1];
			for (j = 0; j < nv; ++j)
				comtr[j] = 0;
			comtr[perm[i]] = 1;
/* transform basis and transformation matrix */
			change(*b, rank, tr, dim);
			change(com->trans, rank, tr, nv);
			*b = (int**)xrealloc(*b, (rank+1) * sizeof(int*));
			(*b)[rank] = (int*)xmalloc(dim * sizeof(int));
			com->trans = (int**)xrealloc(com->trans, (rank+1) * sizeof(int*));
			com->trans[rank] = (int*)xmalloc(nv * sizeof(int));
			f = (int**)xrealloc(f, (rank+1) * sizeof(int*));
			for (j = 0; j < rank; ++j)
			{
				f[j] = (int*)xrealloc(f[j], (rank+1) * sizeof(int));
			}
			f[rank] = (int*)xmalloc((rank+1) * sizeof(int));
			for (j = 0; j < rank; ++j)
			{
				fj = f[j];
				for (k = 0; k <= j; ++k)
				{
					fj[k] = scp((*b)[j], F, (*b)[k], dim);
					f[k][j] = fj[k];
				}
			}
			tr = (int**)xrealloc(tr, (rank+1) * sizeof(int*));
			for (j = 0; j < rank; ++j)
			{
				tr[j] = (int*)xrealloc(tr[j], (rank+1) * sizeof(int));
			}
			tr[rank] = (int*)xmalloc((rank+1) * sizeof(int));
		}
		else if (abs(tr[rank][rank]) > 1)
/* v[perm[i]] enlarges the lattice generated by the vectors in b by a finite
   index |tr[rank][rank]| */
		{
			for (j = 0; j < dim; ++j)
				(*b)[rank][j] = vi[j];
			comtr = com->trans[rank];
			for (j = 0; j < nv; ++j)
				comtr[j] = 0;
			comtr[perm[i]] = 1;
/* transform basis and transformation matrix */
			change(*b, rank+1, tr, dim);
			change(com->trans, rank+1, tr, nv);
			for (j = 0; j < rank; ++j)
			{
				fj = f[j];
				for (k = 0; k <= j; ++k)
				{
					fj[k] = scp((*b)[j], F, (*b)[k], dim);
					f[k][j] = fj[k];
				}
			}
		}
		++i;
	}
	com->rank = rank;
	for (i = 0; i < nv; ++i)
		free(Fv[i]);
	free(Fv);
	for (i = 0; i <= rank; ++i)
	{
		free(f[i]);
		free(tr[i]);
	}
	free(f);
	free(tr);
	free(perm);
}

/************************************************************************\
|	expresses the vectors in v in the basis b,
|	puts the transformation matrix on com->coef, i.e. v = com->coef * b,	
|	uses LLL-reduction with respect to F to obtain the coefficients	
\************************************************************************/
void 
coef (scpcomb *com, int **b, int **v, int **F, int dim)
{
	int	i, j, **Fb, nb, **Fv, nv, **f, **tr, sign;
	int	*fi, *fnb, *trnb, *Fvi, *comci;

	com->coef = (int**)xmalloc((com->list.n+1) * sizeof(int*));
	for (i = 0; i <= com->list.n; ++i)
	{
		com->coef[i] = (int*)xmalloc(com->rank * sizeof(int));
	}
	nb = com->rank;
	nv = com->list.n + 1;
/* Fv is the list of products of F with the vectors in v, scalar products
   with respect to F can be calculated as standard scalar products of
   vectors in v with vectors in Fv */
	Fv = (int**)xmalloc(nv * sizeof(int*));
	for (i = 0; i < nv; ++i)
	{
		Fv[i] = (int*)xmalloc(dim * sizeof(int));
		for (j = 0; j < dim; ++j)
			Fv[i][j] = sscp(F[j], v[i], dim);
	}
/* Fb is the list of products of F with the vectors in b */
	Fb = (int**)xmalloc(nb * sizeof(int*));
	for (i = 0; i < nb; ++i)
	{
		Fb[i] = (int*)xmalloc(dim * sizeof(int));
		for (j = 0; j < dim; ++j)
			Fb[i][j] = sscp(F[j], b[i], dim);
	}
	f = (int**)xmalloc((nb+1) * sizeof(int*));
	for (i = 0; i <= nb; ++i)
	{
		f[i] = (int*)xmalloc((nb+1) * sizeof(int));
	}
	for (i = 0; i < nb; ++i)
	{
		fi = f[i];
		for (j = 0; j <= i; ++j)
		{
			fi[j] = sscp(b[i], Fb[j], dim);
			f[j][i] = fi[j];
		}
	}
	tr = (int**)xmalloc((nb+1) * sizeof(int*));
	for (i = 0; i <= nb; ++i)
	{
		tr[i] = (int*)xmalloc((nb+1) * sizeof(int));
	}
	fnb = f[nb];
	trnb = tr[nb];
	for (i = 0; i < nv; ++i)
	{
		Fvi = Fv[i];
		for (j = 0; j < nb; ++j)
		{
			f[j][nb] = sscp(b[j], Fvi, dim);
			fnb[j] = f[j][nb];
		}
		fnb[nb] = sscp(v[i], Fvi, dim);
/* if the vector v[i] is in the lattice generated by the base b, the rank of the
   Gram-matrix f must be nb and the last row of the transformation matrix tr
   expresses the 0-vector, in particular |tr[nb][nb]| must be <=1, since 
   otherwise the vector v[i] would enlarge the lattice be an index 
   |tr[nb][nb]| */
		if (lll(f, tr, nb+1) > nb  ||  abs(trnb[nb]) > 1)
		{
			fprintf(stderr, "Error: vector lies not in lattice\n");
			exit (3);
		}
		else
		{
			comci = com->coef[i];
			sign = trnb[nb];
			for (j = 0; j < nb; ++j)
				comci[j] = -sign * trnb[j];
		}
	}
	for (i = 0; i < nv; ++i)
		free(Fv[i]);
	free(Fv);
	for (i = 0; i < nb; ++i)
		free(Fb[i]);
	free(Fb);
	for (i = 0; i <= nb; ++i)
	{
		free(f[i]);
		free(tr[i]);
	}
	free(f);
	free(tr);
}

/************************************************************************\
|   com->F[i] is the Gram-matrix of the basis b with respect to F.A[i]
\************************************************************************/
void 
scpforms (scpcomb *com, int **b, invar F)
{
	int	i, j, k, dim, **Fbi, nb, **FAi, *comFij;

	dim = F.dim;
	nb = com->rank;
	com->F = (int***)xmalloc(F.n * sizeof(int**));
/* Fbi is the list of products of F.A[i] with the vectors in b */
	Fbi = (int**)xmalloc(nb * sizeof(int*));
	for (j = 0; j < nb; ++j)
	{
		Fbi[j] = (int*)xmalloc(dim * sizeof(int));
	}
	for (i = 0; i < F.n; ++i)
	{
		com->F[i] = (int**)xmalloc(nb * sizeof(int*));
		FAi = F.A[i];
		for (j = 0; j < nb; ++j)
		{
			for (k = 0; k < dim; ++k)
				Fbi[j][k] = sscp(FAi[k], b[j], dim);
		}
		for (j = 0; j < nb; ++j)
		{
			com->F[i][j] = (int*)xmalloc(nb * sizeof(int));
			comFij = com->F[i][j];
			for (k = 0; k < nb; ++k)
				comFij[k] = sscp(b[j], Fbi[k], dim);
		}
	}
	for (j = 0; j < nb; ++j)
		free(Fbi[j]);
	free(Fbi);
}
