/*****	This file contains some routines for orbit calculations	*****/
#include"typedef.h"

/**********************************************************************\
|	V.v is a sorted list of length V.n of vectors
|	of dimension V.dim, the number of V.v[nr]*A in
|	the list is returned, where a negative number 
|	indicates the negative of a vector
\**********************************************************************/
static int operate(nr, A, V)
veclist	V;
int	nr, **A;
{
	int	i, im, *w;

	if ((w = (int*)malloc(V.dim * sizeof(int))) == 0)
		exit (2);
	vecmatmul(V.v[abs(nr)], A, V.dim, w);
	if (nr < 0)
	{
		for (i = 0; i < V.dim; ++i)
			w[i] *= -1;
	}
	im = numberof(w, V);
	if (abs(im) > V.n)
/* the vector is not in the list */
	{
		fprintf(stderr, "Error: image of vector %d not found\n", nr);
		exit (3);
	}
	free(w);
	return(im);
}

/**********************************************************************\
|	computes the orbit of npt points in pt 
|	under the nG matrices in G and puts the
|	orbit in orb, allocates the memory for 
|	orb, returns the length of the orbit, 
|	the points are the indices in the list V.v
\**********************************************************************/
static int orbit(pt, npt, G, nG, V, orb)
veclist	V;
int	*pt, npt, ***G, nG, **orb;
{
	int	i, norb, cnd, im, *flag;

	if ((*orb = (int*)malloc(npt * sizeof(int))) == 0)
		exit (2);
/* if flag[i + V.n] = 1, then the point i is already in the orbit */
	if ((flag = (int*)malloc((2*V.n + 1) * sizeof(int))) == 0)
		exit (2);
	for (i = 0; i <= 2*V.n; ++i)
		flag[i] = 0;
	for (i = 0; i < npt; ++i)
	{
		(*orb)[i] = pt[i];
		flag[pt[i]+V.n] = 1;
	}
	norb = npt;
	cnd = 0;
	while (cnd < norb)
	{
		for (i = 0; i < nG; ++i)
		{
			im = operate((*orb)[cnd], G[i], V);
			if (flag[im+V.n] == 0)
/* the image is a new point in the orbit */
			{
				++norb;
				if ((*orb = (int*)realloc(*orb, norb * sizeof(int))) == 0)
					exit (2);
				(*orb)[norb-1] = im;
				flag[im+V.n] = 1;
			}
		}
		++cnd;
	}
	free(flag);
	return(norb);
}

/**********************************************************************\
|	checks, whether the orbit of pt under 
|	the nG matrices in G has at least length orblen
\**********************************************************************/
static int orbitlen(pt, orblen, G, nG, V)
veclist	V;
int	pt, orblen, ***G, nG;
{
	int	i, len, cnd, im, *orb, *flag;

	if ((orb = (int*)malloc(orblen * sizeof(int))) == 0)
		exit(2);
/* if flag[i + V.n] = 1, then the point i is already in the orbit */
	if ((flag = (int*)malloc((2*V.n + 1) * sizeof(int))) == 0)
		exit(2);
	for (i = 0; i <= 2*V.n; ++i)
		flag[i] = 0;
	orb[0] = pt;
	flag[pt+V.n] = 1;
	for (i = 1; i < orblen; ++i)
		orb[i] = 0;
	len = 1;
	cnd = 0;
	while (cnd < len  &&  len < orblen)
	{
		for (i = 0; i < nG  &&  len < orblen; ++i)
		{
			im = operate(orb[cnd], G[i], V);
			if (flag[im+V.n] == 0)
/* the image is a new point in the orbit */
			{
				orb[len] = im;
				++len;
				flag[im+V.n] = 1;
			}
		}
		++cnd;
	}
	free(orb);
	free(flag);
	return(len);
}

/**********************************************************************\
|	deletes the elements in orb2 from orb1, 
|	an entry 0 marks the end of the list, returns the length of orb1
\**********************************************************************/
static int delete(orb1, l1, orb2, l2)
int	*orb1, l1, *orb2, l2;
{
	int	i, j, len, o2i;

/* the true length of orb1 is len, hence orb1[len-1] != 0 */
	for (i = 0; i < l1  &&  orb1[i] != 0; ++i);
	len = i;
	for (i = 0; i < l2  &&  orb2[i] != 0; ++i)
	{
		o2i = orb2[i];
		for (j = 0; j < len  &&  orb1[j] != o2i; ++j);
/* orb1[j] = orb2[i], hence delete orb1[j] from orb1 */
		if (j < len)
		{
			orb1[j] = orb1[len-1];
			orb1[len-1] = 0;
			--len;
		}	
	}
	return(len);
}

/**********************************************************************\
|	computes the orbit of fp.e[I] under the 
|	generators in G->g[I]...G->g[n-1] and elements 
|	stabilizing fp.e[I],
|	has some heuristic break conditions,
|	the generators in G->g[i] stabilize 
|	fp.e[0]...fp.e[i-1] but not fp.e[i], 
|	G->ng[i] is the number of generators in G->g[i],
|	the first G->nsg[i] of which are elements which
|	are obtained as stabilizer elements in 
|	<G->g[0],...,G->g[i-1]>, G->ord[i] is the orbit
|	length of fp.e[i] under <G->g[i],...,G->g[n-1]>
\**********************************************************************/
static void stab(I, G, fp, V)
int I;
group	*G;
fpstruct fp;
veclist	V;
{
	int	*orb, len, cnd, tmplen;
	int	**w, *flag, ***H, ***Hj, **S, **tmp, ***Ggj;
	int	i, j, k, l, dim, im, nH, nHj, fail;
	int	Maxfail, Rest;

/* some heuristic break conditions for the computation of stabilizer elements:
   it would be too expensive to calculate all the stabilizer generators, which
   are obtained from the orbit, since this is highly redundant, 
   on the other hand every new generator which enlarges the group is much 
   cheaper than one obtained from the backtrack,
   after Maxfail subsequent stabilizer elements, that do not enlarge the group,
   Rest more elements are calculated even if they leave the group unchanged,
   since it turned out that this is often useful in the following steps,
   increasing the parameters will possibly decrease the number of generators
   for the group, but will increase the running time,
   there is no magic behind this heuristic, tuning might be appropriate */
	dim = V.dim;
	Rest = 0;
	for (i = I; i < dim; ++i)
	{
		if (fp.diag[i] > 1  &&  G->ord[i] < fp.diag[i])
			++Rest;
	}
	Maxfail = Rest;
	for (i = 0; i < dim; ++i)
	{
		if (fp.diag[i] > 1)
			++Maxfail;
	}
	nH = 0;
	for (i = I; i < dim; ++i)
		nH += G->ng[i];
/* H are the generators of the group in which the stabilizer is computed */
	if ((H = (int***)malloc(nH * sizeof(int**))) == 0)
		exit (2);
	if ((Hj = (int***)malloc((nH+1) * sizeof(int**))) == 0)
		exit (2);
	k = 0;
	for (i = I; i < dim; ++i)
	{
		for (j = 0; j < G->ng[i]; ++j)
		{
			H[k] = G->g[i][j];
			++k;
		}
	}
/* in w[V.n+i] an element is stored that maps fp.e[I] on v[i] */
	if ((w = (int**)malloc((2*V.n+1) * sizeof(int*))) == 0)
		exit (2);
/* orb contains the orbit of fp.e[I] */
	if ((orb = (int*)malloc(2*V.n * sizeof(int))) == 0)
		exit (2);
	for (i = 0; i < 2*V.n; ++i)
		orb[i] = 0;
/* if flag[i + V.n] = 1, then the point i is already in the orbit */
	if ((flag = (int*)malloc((2*V.n+1) * sizeof(int))) == 0)
		exit (2);
	for (i = 0; i <= 2*V.n; ++i)
		flag[i] = 0;
/* S is a matrix to hold a stabilizer element temporarily */
	if ((S = (int**)malloc(dim * sizeof(int*))) == 0)
		exit (2);
	for (i = 0; i < dim; ++i)
	{
		if ((S[i] = (int*)malloc(dim * sizeof(int))) == 0)
			exit (2);
	}
	orb[0] = fp.e[I];
	flag[orb[0]+V.n] = 1;
	if ((w[orb[0]+V.n] = (int*)malloc(dim * sizeof(int))) == 0)
		exit (2);
	for (i = 0; i < dim; ++i)
		w[orb[0]+V.n][i] = fp.e[i];
	cnd = 0;
	len = 1;
/* fail is the number of successive failures */
	fail = 0;
	while (len > cnd  &&  fail < Maxfail+Rest)
	{
		for (i = 0; i < nH  &&  fail < Maxfail+Rest; ++i)
		{
			if (fail >= Maxfail)
/* there have already been Maxfail successive failures, now a random generator
   is applied to a random point of the orbit to get Rest more stabilizer 
   elements */
			{
				cnd = rand() % len;
				if (cnd < 0)
					cnd += len;
				i = rand() % nH;
				if (i < 0)
					i += nH;
			}
			im = operate(orb[cnd], H[i], V);
			if (flag[im+V.n] == 0)
/* a new element is found, appended to the orbit and an element mapping
   fp.e[I] to im is stored in w[im+V.n] */
			{
				orb[len] = im;
				flag[im+V.n] = 1;
				if ((w[im+V.n] = (int*)malloc(dim * sizeof(int))) == 0)
					exit (2);
				for (j = 0; j < dim; ++j)
					w[im+V.n][j] = operate(w[orb[cnd]+V.n][j], H[i], V);
				++len;
			}
			else
/* the image was already in the orbit */
			{
/* j is the first index where the images of the old and the new element
   mapping e[I] on im differ */
				for (j = I; j < dim  &&  operate(w[orb[cnd]+V.n][j], H[i], V) == w[im+V.n][j]; ++j);
				if (j < dim  &&  
				    (G->ord[j] < fp.diag[j] || fail >= Maxfail))
				{
/* new stabilizer element S = w[orb[cnd]+V.n] * H[i] * (w[im+V.n])^-1 */
					stabil(S, w[orb[cnd]+V.n], w[im+V.n], fp.per, H[i], V);
					Hj[0] = S;
					nHj = 1;
					for (k = j; k < dim; ++k)
					{
						for (l = 0; l < G->ng[k]; ++l)
						{
							Hj[nHj] = G->g[k][l];
							++nHj;
						}
					}
					tmplen = orbitlen(fp.e[j], fp.diag[j], Hj, nHj, V);
					if (tmplen > G->ord[j]  ||  fail >= Maxfail)
/* the new stabilizer element S enlarges the orbit of e[j] */
					{
						++G->ng[j];
						++G->nsg[j];
/* allocate memory for the new generator */
						if ((G->g[j] = (int***)realloc(G->g[j], G->ng[j] * sizeof(int**))) == 0)
							exit (2);
						if ((G->g[j][G->ng[j]-1] = (int**)malloc(dim * sizeof(int*))) == 0)
							exit (2);
						for (k = 0; k < dim; ++k)
						{
							if ((G->g[j][G->ng[j]-1][k] = (int*)malloc(dim * sizeof(int))) == 0)
								exit (2);
						}
						Ggj = G->g[j];
/* the new generator is inserted as stabilizer element nr. nsg[j]-1 */
						for (k = G->ng[j]-1; k >= G->nsg[j]; --k)
						{
							tmp = Ggj[k];
							Ggj[k] = Ggj[k-1];
							Ggj[k-1] = tmp;
						}
						for (k = 0; k < dim; ++k)
						{
							for (l = 0; l < dim; ++l)
								Ggj[G->nsg[j]-1][k][l] = S[k][l];
						}
						G->ord[j] = tmplen;
						++nH;
						if ((H = (int***)realloc(H, nH * sizeof(int**))) == 0)
							exit (2);
						if ((Hj = (int***)realloc(Hj, (nH+1) * sizeof(int**))) == 0)
							exit (2);
/* the new generator is appended to H */
						H[nH-1] = Ggj[G->nsg[j]-1];
/* the number of failures is reset to 0 */ 
						if (fail < Maxfail)
							fail = 0;
						else
							++fail;
					}
					else
/* the new stabilizer element S does not enlarge the orbit of e[j] */
						++fail;
				}
				else if (j < dim  &&  fail < Maxfail  ||  
					 j == dim  &&  fail >= Maxfail)
					++fail;
/* if S is the identity and fail < Maxfail, nothing is done */
			}
		}
		if (fail < Maxfail)
			++cnd;
	}
	for (i = 0; i <= 2*V.n; ++i)
	{
		if (flag[i] == 1)
			free(w[i]);
	}
	free(w);
	for (i = 0; i < dim; ++i)
		free(S[i]);
	free(S);
	free(orb);
	free(flag);
	free(H);
	free(Hj);
}

/**********************************************************************\
|	generates the matrix X which has as row
|	per[i] the vector nr. x[i] from the list v
\**********************************************************************/
static void matgen(x, X, dim, per, v)
int	*x, **X, dim, *per, **v;
{
	int	i, j, xi, *Xperi;

	for (i = 0; i < dim; ++i)
	{
		Xperi = X[per[i]];
		if ((xi = x[i]) > 0)
		{
			for (j = 0; j < dim; ++j)
				Xperi[j] = v[xi][j];
		}
		else
		{
			for (j = 0; j < dim; ++j)
				Xperi[j] = -v[-xi][j];
		}
	}
}

/**********************************************************************\
|	x1 corresponds to an element X1 mapping some vector e on p1, 
|	x2 to an element X2 mapping e on p2 and G is a generator mapping 
|	p1 on p2, then S = X1*G*X2^-1 stabilizes e
\**********************************************************************/
static void stabil(S, x1, x2, per, G, V)
veclist	V;
int	**S, *x1, *x2, *per, **G;
{
	int	i, dim, *x, **XG, **X2;

	dim = V.dim;
	if ((XG = (int**)malloc(dim * sizeof(int*))) == 0)
		exit (2);
	if ((X2 = (int**)malloc(dim * sizeof(int*))) == 0)
		exit (2);
	for (i = 0; i < dim; ++i)
	{
		if ((XG[i] = (int*)malloc(dim * sizeof(int))) == 0)
			exit (2);
		if ((X2[i] = (int*)malloc(dim * sizeof(int))) == 0)
			exit (2);
	}
	if ((x = (int*)malloc(dim * sizeof(int))) == 0)
		exit (2);
	for (i = 0; i < dim; ++i)
		x[i] = operate(x1[i], G, V);
/* XG is the composite mapping of the matrix corresponding to x1 and G */
	matgen(x, XG, dim, per, V.v);
	matgen(x2, X2, dim, per, V.v);
/* S = XG * X2^-1 */
	psolve(S, X2, XG, dim, V.prime);
	free(x);
	for (i = 0; i < dim; ++i)
	{
		free(XG[i]);
		free(X2[i]);
	}
	free(XG);
	free(X2);
}
