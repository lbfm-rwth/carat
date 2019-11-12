/*****	This file contains some routines to compute automorphisms of a lattice 
	and to determine candidates for the image of a base vector	*****/
#include "typedef.h"
#include "utils.h"
#include "types.h"


/***************************************************************\
|   tests, whether x[0],...,x[I-1] is a partial 
|   automorphism, using scalar product combinations
|   and Bacher-polynomials depending on the chosen 
|   options, puts the candidates for x[I] (i.e. the
|   vectors vec such that the scalar products of 
|   x[0],...,x[I-1],vec are correct) on CI, 
|   returns their number (should be fp[I])
\***************************************************************/
int 
cand (int *CI, int I, int *x, veclist V, invar F, fpstruct fp, scpcomb *comb, bachpol *bach, flagstruct flags)
{
	int	i, j, k, dim, okp, okm, sign, nr, fail, num;
	int	*vec, *scpvec, **xvec, **xbase, **Fxbase, DEP, len, rank, n;
	int	*Vvj, *FAiI, **Fvi, *xnum, xk, *comtri, *comcoi, xbij, vj;
	scpcomb	com;
        int test;

	dim = F.dim;
	DEP = flags.DEPTH;
	len = F.n * DEP;

        if(normal_option == TRUE && I >0)
        {
           test = normal_aut_test(x, I-1, V);
           if(test == FALSE)
             return(0);
        }
	if (I-1 >= 0  &&  I-1 < flags.BACHDEP)
	{
		if (flags.BACH[2] == 0)
			flags.BACHSCP = V.v[abs(x[I-1])][dim] / 2;
		if (bachcomp(bach[I-1], x[I-1], flags.BACHSCP, V, F.v[0]) == 0)
			return(0);
	}
	vec = (int*)xmalloc(dim * sizeof(int));
	scpvec = (int*)xmalloc(len * sizeof(int));
	if (I-1 >= 0  &&  DEP > 0)
	{
		com = comb[I-1];
		rank = com.rank;
		n = com.list.n;
/* xvec is the list of vector sums which are computed with respect to the 
   partial basis in x */
		xvec = (int**)xmalloc((n+1) * sizeof(int*));
		for (i = 0; i <= n; ++i)
		{
			xvec[i] = (int*)xmalloc(dim * sizeof(int));
			for (j = 0; j < dim; ++j)
				xvec[i][j] = 0;
		}
/* xbase should be a basis for the lattice generated by the vectors in xvec,
   it is obtained via the transformation matrix comb[I-1].trans */
		xbase = (int**)xmalloc(rank * sizeof(int*));
		for (i = 0; i < rank; ++i)
		{
			xbase[i] = (int*)xmalloc(dim * sizeof(int));
		}
/* Fxbase is the product of a form F with the base xbase */
		Fxbase = (int**)xmalloc(rank * sizeof(int*));
		for (i = 0; i < rank; ++i)
		{
			Fxbase[i] = (int*)xmalloc(dim * sizeof(int));
		}
	}
	else
	{
	    memset(&com, 0, sizeof(com));
	    xvec = 0;
	    xbase = 0;
	    Fxbase = 0;
	    rank = 0;
	    n = 0;
	}
/* CI is the list for the candidates */
	for (i = 0; i < fp.diag[I]; ++i)
		CI[i] = 0;
	nr = 0;
	fail = 0;
	for (j = 1; j <= V.n  &&  fail == 0; ++j)
	{
		Vvj = V.v[j];
		okp = 0;
		okm = 0;
		for (k = 0; k < len; ++k)
			scpvec[k] = 0;
		for (i = 0; i < F.n; ++i)
		{
			FAiI = F.A[i][fp.per[I]];
			Fvi = F.v[i];
/* vec is the vector of scalar products of V.v[j] with the first I base vectors
   x[0]...x[I-1] */
			for (k = 0; k < I; ++k)
			{
				if ((xk = x[k]) > 0)
					vec[k] = sscp(Vvj, Fvi[xk], dim);
				else
					vec[k] = -sscp(Vvj, Fvi[-xk], dim);
			}
			for (k = 0; k < I  &&  vec[k] == FAiI[fp.per[k]]; ++k);
			if (k == I  &&  Vvj[dim+i] == FAiI[fp.per[I]])
/* V.v[j] is a candidate for x[I] with respect to the form F.A[i] */
				++okp;
			for (k = 0; k < I  &&  vec[k] == -FAiI[fp.per[k]]; ++k);
			if (k == I  &&  Vvj[dim+i] == FAiI[fp.per[I]])
/* -V.v[j] is a candidate for x[I] with respect to the form F.A[i] */
				++okm;
			if (I-1 >= 0  &&  DEP > 0)
			{
				for (k = I-1; k >= 0  &&  k > I-1-DEP; --k)
					scpvec[i*DEP+I-1-k] = vec[k];
			}
		}
		if (I-1 >= 0  &&  DEP > 0)
/* check, whether the scalar product combination scpvec is contained in the
   list comb[I-1].list */
		{
			for (i = 0; i < len  &&  scpvec[i] == 0; ++i);
			if (i == len)
				num = 0;
			else
			{
				num = numberof(scpvec, com.list);
				sign = num > 0 ? 1 : -1;
				num *= sign;
			}
			if (num > n)
/* scpvec is not found, hence x[0]...x[I-1] is not a partial automorphism */
				fail = 1;
			else if (num > 0)
/* scpvec is found and the vector is added to the corresponding vector sum */
			{
				xnum = xvec[num];
				for (k = 0; k < dim; ++k)
					xnum[k] += sign * Vvj[k];
			}
		}
		if (okp == F.n)
/* V.v[j] is a candidate for x[I] */
		{
			if (nr < fp.diag[I])
			{
				CI[nr] = j;
				++nr;
			}
			else
/* there are too many candidates */
				fail = 1;
		}
		if (okm == F.n)
/* -V.v[j] is a candidate for x[I] */
		{
			if (nr < fp.diag[I])
			{
				CI[nr] = -j;
				++nr;
			}
			else
/* there are too many candidates */
				fail = 1;
		}
	}
	if (fail == 1)
		nr = 0;
	if (nr == fp.diag[I]  &&  I-1 >= 0  &&  DEP > 0)
/* compute the basis of the lattice generated by the vectors in xvec via the
   transformation matrix comb[I-1].trans */
	{
		for (i = 0; i < rank; ++i)
		{
			comtri = com.trans[i];
			for (j = 0; j < dim; ++j)
			{
				xbij = 0;
				for (k = 0; k <= n; ++k)
					xbij += comtri[k] * xvec[k][j];
				xbase[i][j] = xbij;
			}
		}
	}
	if (nr == fp.diag[I]  &&  I-1 >= 0  &&  DEP > 0)
/* check, whether the base xbase has the right scalar products */
	{
		for (i = 0; i < F.n  &&  nr > 0; ++i)
		{
			for (j = 0; j < rank; ++j)
			{
				for (k = 0; k < dim; ++k)
					Fxbase[j][k] = sscp(F.A[i][k], xbase[j], dim);
			}
			for (j = 0; j < rank  &&  nr > 0; ++j)
			{
				for (k = 0; k <= j  &&  nr > 0; ++k)
				{
					if (sscp(xbase[j], Fxbase[k], dim) != com.F[i][j][k])
/* a scalar product is wrong */
						nr = 0;
				}
			}
		}
	}
	if (nr == fp.diag[I]  &&  I-1 >= 0  &&  DEP > 0)
/* check, whether comb[I-1].coef * xbase = xvec */
	{
		for (i = 0; i <= n  &&  nr > 0; ++i)
		{
			comcoi = com.coef[i];
			for (j = 0; j < dim; ++j)
			{
				vj = 0;
				for (k = 0; k < rank; ++k)
					vj += comcoi[k] * xbase[k][j];
				if (vj != xvec[i][j])
/* an entry is wrong */
					nr = 0;
			}
		}
	}
	if (I-1 >= 0  &&  DEP > 0)
	{
		for (i = 0; i <= n; ++i)
			free(xvec[i]);
		free(xvec);
		for (i = 0; i < rank; ++i)
		{
			free(xbase[i]);
			free(Fxbase[i]);
		}
		free(xbase);
		free(Fxbase);
	}
	free(scpvec);
	free(vec);
	return(nr);
}

/**************************************************\
|   search new automorphisms until on all
|   levels representatives for all orbits
|   have been tested
\**************************************************/
void 
autom (group *G, veclist V, invar F, fpstruct fp, scpcomb *comb, bachpol *bach, flagstruct flags)
{
	FILE	*outfile;
	int	i, j, dim, step, im, **C, nC, ***H, nH, **Ggng, found, tries;
	int	*x, *orb, *oc, noc, *bad, nbad;

	dim = F.dim;
	C = (int**)xmalloc(dim * sizeof(int*));
/* C[i] is the list of candidates for the image of the i-th base-vector */
	for (i = 0; i < dim; ++i)
	{
		C[i] = (int*)xmalloc(fp.diag[i] * sizeof(int));
	}
/* x is the new base, i.e. x[i] is the index in V.v of the i-th base-vector */
	x = (int*)xmalloc(dim * sizeof(int));
	bad = (int*)xmalloc(2*V.n * sizeof(int));
	for (step = flags.STAB; step < dim; ++step)
	{
		nH = 0;
		for (i = step; i < dim; ++i)
			nH += G->ng[i];
		H = (int***)xmalloc(nH * sizeof(int**));
		nH = 0;
		for (i = step; i < dim; ++i)
		{
			for (j = 0; j < G->ng[i]; ++j)
			{
				H[nH] = G->g[i][j];
				++nH;
			}
		}
		for (i = 0; i < 2*V.n; ++i)
			bad[i] = 0;
		nbad = 0;
/* the first step base-vectors are fixed */
		for (i = 0; i < step; ++i)
			x[i] = fp.e[i];
/* compute the candidates for x[step] */
		if (fp.diag[step] > 1)
/* if fp.diag[step] compute the candidates for x[step] */
			nC = cand(C[step], step, x, V, F, fp, comb, bach, flags);
		else
/* if fp.diag[step] == 1, fp.e[step] is the only candidate */
		{
			C[step][0] = fp.e[step];
			nC = 1;
		}
		G->ord[step] = orbit(&(fp.e[step]), 1, H, nH, V, &orb);
/* delete the orbit of the step-th base-vector from the candidates */
		nC = delete_from_orbit(C[step], nC, orb, G->ord[step]);
		free(orb);
		while (nC > 0  &&  (im = C[step][0]) != 0)
		{
/* try vector V.v[im] as image of the step-th base-vector */
			x[step] = im;
			if (step < dim-1)
			{
/* check, whether x[0]...x[step] is a partial basis and compute the candidates
   for x[step+1] */
				if (cand(C[step+1], step+1, x, V, F, fp, comb, bach, flags) == fp.diag[step+1])
/* go into the recursion */
					found = aut(step+1, x, C, G, V, F, fp, comb, bach, flags);
				else
					found = 0;
			}
			else
                        {
                           found = 1;
                           if(normal_option == TRUE)
                           {
                           found = 0;
                           for(i=0;i<fp.diag[step] && C[step][0] != 0 && found == 0; i++)
                           {
                              x[step] = C[step][0];
                              found = normal_aut_test(x, step, V);
                              if(found == 0)
                              {
                                 for(j=0;j<fp.diag[step]-1;j++)
                                  C[step][j] = C[step][j+1];
                                 C[step][fp.diag[step]-1] = 0;
                              }
                           }
                           }
                        }
			if (found == 0)
/* x[0]...x[step] can not be continued to an automorphism */
			{
				noc = orbit(&im, 1, H, nH, V, &oc);
/* delete the orbit of im from the candidates for x[step] */
				nC = delete_from_orbit(C[step], nC, oc, noc);
				free(oc);
				bad[nbad] = im;
				++nbad;
			}
			else
/* a new generator has been found */
			{
				++G->ng[step];
				G->g[step] = (int***)xrealloc(G->g[step], G->ng[step] * sizeof(int**));
				G->g[step][G->ng[step]-1] = (int**)xmalloc(dim * sizeof(int*));
				for (i = 0; i < dim; ++i)
				{
					G->g[step][G->ng[step]-1][i] = (int*)xmalloc(dim * sizeof(int));
				}
				Ggng = G->g[step][G->ng[step]-1];
/* append the new generator to G->g[step] */
				matgen(x, Ggng, dim, fp.per, V.v);
				++nH;
				H = (int***)xrealloc(H, nH * sizeof(int**));
				nH = 0;
				for (i = step; i < dim; ++i)
				{
					for (j = 0; j < G->ng[i]; ++j)
					{
						H[nH] = G->g[i][j];
						++nH;
					}
				}
				if (flags.PRINT == 1)
/* if the -P option is given, print the new generator on AUTO.tmp */
				{
					outfile = fopen("AUTO.tmp", "a");
					fprintf(outfile, "%d\n", dim);
					for (i = 0; i < dim; ++i)
					{
						for (j = 0; j < dim; ++j)
							fprintf(outfile, " %2d", Ggng[i][j]);
						fprintf(outfile, "\n");
					}
					fclose(outfile);
				}
/* compute the new orbit of fp.e[step] */
				G->ord[step] = orbit(&(fp.e[step]), 1, H, nH, V, &orb);
/* delete the orbit from the candidates for x[step] */
				nC = delete_from_orbit(C[step], nC, orb, G->ord[step]);
				free(orb);
/* compute the new orbit of the vectors, which could not be continued to an
   automorphism */
				noc = orbit(bad, nbad, H, nH, V, &oc);
/* delete the orbit from the candidates */
				nC = delete_from_orbit(C[step], nC, oc, noc);
				free(oc);
			}
		}
/* test, whether on step flags.STAB some generators may be omitted */
		if (step == flags.STAB)
		{
			for (tries = G->nsg[step]; tries < G->ng[step]; ++tries)
			{
				nH = 0; 
				for (j = 0; j < tries; ++j)
				{
					H[nH] = G->g[step][j];
					++nH;
				}
				for (j = tries+1; j < G->ng[step]; ++j)
				{
					H[nH] = G->g[step][j];
					++nH;
				}
				for (i = step+1; i < dim; ++i)
				{
					for (j = 0; j < G->ng[i]; ++j)
					{
						H[nH] = G->g[i][j];
						++nH;
					}
				}
				if (orbitlen(fp.e[step], G->ord[step], H, nH, V) == G->ord[step])
	/* the generator g[step][tries] can be omitted */
				{
					for (i = 0; i < dim; ++i)
						free(G->g[step][tries][i]);
					free(G->g[step][tries]);
					for (i = tries; i < G->ng[step]-1; ++i)
						G->g[step][i] = G->g[step][i+1];
					--G->ng[step];
					--tries;
				}
			}
		}
		free(H);
		if (step < dim-1  &&  G->ord[step] > 1)
/* calculate stabilizer elements fixing the basis-vectors 
   fp.e[0]...fp.e[step] */
			stab(step, G, fp, V);
		if (flags.PRINT == 1)
/* if the -P option is given, print G.ord[step] AUTO.tmp */
		{
			outfile = fopen("AUTO.tmp", "a");
			fprintf(outfile, "G.ord[%d] = %d\n", step, G->ord[step]);
			fclose(outfile);
		}
	}
	for (i = 0; i < dim; ++i)
		free(C[i]);
	free(C);
	free(x);
	free(bad);
}

/***************************************************************\
|	the heart of the program: the recursion	
\***************************************************************/
int 
aut (int step, int *x, int **C, group *G, veclist V, invar F, fpstruct fp, scpcomb *comb, bachpol *bach, flagstruct flags)
{
	int	i, j, dim, orb[1], found;

	dim = F.dim;
	found = 0;
	while (C[step][0] != 0  &&  found == 0)
	{
		if (step < dim-1)
		{
/* choose the image of the base-vector nr. step */
			x[step] = C[step][0];
/* check, whether x[0]...x[step] is a partial automorphism and compute the
   candidates for x[step+1] */
			if (cand(C[step+1], step+1, x, V, F, fp, comb, bach, flags) == fp.diag[step+1])
/* go deeper into the recursion */
				found = aut(step+1, x, C, G, V, F, fp, comb, bach, flags);
			if (found == 1)
				return(found);
/* delete the chosen vector from the list of candidates */
			orb[0] = x[step];
			delete_from_orbit(C[step], fp.diag[step], orb, 1);
		}
		else
		{
                    if(normal_option == TRUE)
                    {
                           found = 0;
                           for(i=0;i<fp.diag[step] && C[step][0] != 0 && found == 0;i++)
                           {
                              x[step] = C[step][0];
                              found = normal_aut_test(x, step, V);
                              if(found == 0)
                              {
                                 for(j=0;j<fp.diag[step]-1;j++)
                                  C[step][j] = C[step][j+1];
                                 C[step][fp.diag[step]-1] = 0;
                              }
                           }
                    }
		    else
		    {
                        /* a new automorphism is found */
			    x[dim-1] = C[dim-1][0];
			    found = 1;
		    }
		}
	}
	return(found);
}
