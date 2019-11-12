#error THIS FILE IS NOT USED

/*****	Main program for the automorphism program AUTO	*****/

#include "typedef.h"
#include "utils.h"
#include "types.h"
#include "matrix.h"


/*************************************************************************\
| The functions 'autgrp' calculates generators and the order of the group
| G := {g in GL_n(Z) | g * Fo[i] * g^{tr} = Fo[i], 1<= i<= Foanz}
| returned via a pointer to 'bravais_TYP'.
|
| The arguments of autgrp are:
| matrix_TYP	**Fo:		a set of n times n matrices,
|				the first must be positiv definite
| int	 	Foanz:		the number of the matrices given in 'Fo'.
| matrix_TYP	**Erz:		if already element of G are known,
|				they can be used for calculating generators
|				for the whole group.
|				The matrices of known elements can be given 
|				to the function by the pointer 'Erz'.
| int		Erzanz:		The number of matrices given in 'Erz"
| matrix_TYP	*SV:		The rows of the matrix 'SV' must be the vectors
|				x in Z^n with x * F[0] * x^{tr} <= m, where
|				m is the maximal diagonal entry of the Matrix
|				F[0]
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


bravais_TYP *
autgrp (matrix_TYP **Fo, int Foanz, matrix_TYP *SV, matrix_TYP **Erz, int Erzanz, int *options)
{
	FILE		*outfile;
	bachpol		*bach;
	flagstruct	flags;
	scpcomb		*comb;
	group		G;
	invar		F;
	veclist		V, norm;
	fpstruct	fp;
	int		dim, max, fail, *vec;
	int		***H, nH, ngen, **sumveclist, **sumvecbase;
	int		i, j, k, l, n, *Vvi, *Vvj, **Fvi, *Fvij, **FAi;
        bravais_TYP *B;

   normal_option = FALSE;
/* get the flags from the command line */
	getflags(&flags, options);
   if(Erzanz > 0)
     flags.GEN = 1;
   else
     flags.GEN = 0;
/* F.n is the number of invariant forms */
	F.n = Foanz;
	F.A = (int***)xmalloc(F.n * sizeof(int**));
/* read the invariant forms */
        F.dim = Fo[0]->cols;
	dim = F.dim;
        for(i=0;i<F.n;i++)
        {
           if(Fo[i]->cols != F.dim || Fo[i]->rows != F.dim)
           {
              printf("forms in autgrp have different dimension\n");
              exit(2);
           }
        }
        for(i=0;i<F.n;i++)
        {
	  F.A[i] = (int**)xmalloc(dim * sizeof(int*));
          for(j=0;j<F.dim;j++)
          {
	    F.A[i][j] = (int*)xmalloc(dim * sizeof(int));
            for(k=0;k<F.dim;k++)
              F.A[i][j][k] = Fo[i]->array.SZ[j][k];
          }
        }
	if (flags.GEN == 1)
/* some automorphisms are already known */
	{
		ngen = Erzanz;
		H = (int***)xmalloc(ngen * sizeof(int**));
		nH = 0;
		fail = 0;
                for(i=0;i<ngen;i++)
                {
                  if(Erz[i]->cols != dim || Erz[i]->rows != dim)
                  {
                    printf("Elements 'Erz' have wrong dimension\n");
                    exit(3);
                  }
                }
                i = 0;
		while (nH+fail < ngen)
		{
		    H[nH] = (int**)xmalloc(dim * sizeof(int*));
                    for(j=0;j<dim;j++)
                    {
		      H[nH][j] = (int*)xmalloc(dim * sizeof(int));
                      for(k=0;k<dim;k++)
                        H[nH][j][k] = Erz[i]->array.SZ[j][k];
                    }
/* check whether the matrix is really an automorphism, i.e. fixes the forms */
			if (checkgen(H[nH], F) == 0)
/* the matrix is not an automorphism */
			{
				++fail;
				for (j = 0; j < dim; ++j)
					free(H[nH][j]);
                                free(H[nH]);
			}
			else
/* the matrix fixes all forms in F */
				++nH;
                    i++;
		}
	}
	else
		nH = 0;



	V.dim = dim;
	V.len = dim + F.n;
/* get the short vectors */
        V.n = SV->rows;
        V.v = (int**)xmalloc((V.n+1) * sizeof(int*));
        for(i=0;i<=V.n;i++)
        {
          V.v[i] = (int*)xmalloc(V.len * sizeof(int));
        }
        for(i=1;i<= V.n;i++)
        {
            for(j=0;j<=dim;j++)
               V.v[i][j] = SV->array.SZ[i-1][j];
        }
/* the first nonzero entry of each vector is made positive and the maximal
   entry of the vectors is determined */
	max = 1;
	for (i = 1; i <= V.n; ++i)
	{
		Vvi = V.v[i];
		for (j = 0; j < dim  &&  Vvi[j] == 0; ++j);
		if (j < dim  &&  Vvi[j] < 0)
		{
			for (k = j; k < dim; ++k)
				Vvi[k] *= -1;
		}
		for (k = j; k < dim; ++k)
		{
			if (abs(Vvi[k]) > max)
				max = abs(Vvi[k]);
		}
	}
	if (max > MAXENTRY)
/* there is an entry, which is bigger than MAXENTRY, which might cause overflow,
   since the program doesn't use long integer arithmetic */
	{
		fprintf(stderr, "Error: Found entry %d in short vectors.\nTo avoid overflow, the entries should not exceed %d.\n", max, MAXENTRY);
		exit (3);
	}
/* V.prime is a prime p, such that every entry x in the short vectors remains
   unchanged under reduction mod p in symmetric form, i.e. -p/2 < x <= p/2 */
	for (V.prime = 2*max + 1; isprime(V.prime) == 0; ++V.prime);
/* sort the vectors and delete doublets */
	sortvecs(&V);
/* the norm-vector (i.e. the vector of the norms with respect to the different
   invariant forms) of each vector must be equal to the norm-vector of one
   of the vectors from the standard-base */
	norm.v = (int**)xmalloc((dim+1) * sizeof(int*));
	norm.n = dim;
	norm.dim = F.n;
	norm.len = F.n;
	for (i = 1; i <= norm.n; ++i)
	{
/* norm.v[i] is the norm combination of the (i-1)-th base-vector */
		norm.v[i] = (int*)xmalloc(norm.dim * sizeof(int));
		for (k = 0; k < norm.dim; ++k)
			norm.v[i][k] = F.A[k][i-1][i-1];
	}
	sortvecs(&norm);
/* delete those vectors, which can not be the image of any of the standard-base
   vectors */
	checkvecs(&V, F, norm);
	for (i = 1; i <= norm.n; ++i)
		free(norm.v[i]);
	free(norm.v);
/* print the invariant forms, if output is in ASCII-format */
/*********************************************
		printf("#%d\n", F.n);
		for (i = 0; i < F.n; ++i)
			putform(F.A[i], dim);
**********************************************/
/* F.v[i][j] is the transposed of the product of the Gram-matrix F.A[i] with the
   transposed of the vector v[j], hence the scalar product of v[j] and v[k] with
   respect to the Gram-matrix F.A[i] can be computed as standard scalar product
   of v[j] and F.v[i][k] */
	F.v = (int***)xmalloc(F.n * sizeof(int**));
/* the product of the maximal entry in the short vectors with the maximal entry
   in F.v[i] should not exceed MAXNORM to avoid overflow */
	max = MAXNORM / max;
	for (i = 0; i < F.n; ++i)
	{
		FAi = F.A[i];
		F.v[i] = (int**)xmalloc((V.n+1) * sizeof(int*));
		Fvi = F.v[i];
		for (j = 1; j <= V.n; ++j)
		{
			Vvj = V.v[j];
			Fvi[j] = (int*)xmalloc(dim * sizeof(int));
			Fvij = Fvi[j];
			for (k = 0; k < dim; ++k)
			{
				Fvij[k] = sscp(FAi[k], Vvj, dim);
				if (abs(Fvij[k]) > max)
/* some entry in F.v[i] is too large */
				{
					fprintf(stderr, "Error: Found entry %d in F.v.\nTo avoid overflow, the entries should not exceed %d.\n", Fvij[k], max);
					exit (3);
				}
			}
		}
	}
/* fp.per is the order in which the images of the base-vectors are set */
	fp.per = (int*)xmalloc(dim * sizeof(int));
/* fp.e[i] is the index in V.v of the i-th vector of the standard-base in the 
   new order */
	fp.e = (int*)xmalloc(dim * sizeof(int));
/* fp.diag is the diagonal of the fingerprint in the new order, fp.diag[i] is
   an upper bound for the length of the orbit of the i-th base-vector under
   the stabilizer of the preceding ones */
	fp.diag = (int*)xmalloc(dim * sizeof(int));
/* compute the fingerprint */
	fingerprint(&fp, F, V);
	if (flags.PRINT == 1)
/* if the -P option is given, print the diagonal fp.diag and the number of short
   vectors on AUTO.tmp */
	{
		outfile = fopen("AUTO.tmp", "a");
		fprintf(outfile, "%d short vectors\n", 2*V.n);
		fprintf(outfile, "fingerprint diagonal:\n");
		for (i = 0; i < dim; ++i)
			fprintf(outfile, " %2d", fp.diag[i]);
		fprintf(outfile, "\n");
		fclose(outfile);
	}
/* get the standard basis in the new order */
	vec = (int*)xmalloc(dim * sizeof(int));
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
			vec[j] = 0;
		vec[fp.per[i]] = 1;
		fp.e[i] = numberof(vec, V);
		if (abs(fp.e[i]) > V.n)
/* the standard-base must occur in the set of short vectors, since Id is
   certainly an automorphism */
		{
			fprintf(stderr, "Error: standard basis vector nr. %d not found\n", i);
			exit (3);
		}
	}
	free(vec);
/* if the -D option is given, the scalar product combinations and the 
   corresponding vector sums are computed for the standard-basis */
	if (flags.DEPTH > 0)
	{
		comb = (scpcomb*)xmalloc(dim * sizeof(scpcomb));
		for (i = 0; i < dim; ++i)
		{
/* compute the list of scalar product combinations and the corresponding
   vector sums */
			scpvecs(&comb[i].list, &sumveclist, i, fp.e, flags.DEPTH, V, F);
/* compute a basis for the lattice that is generated by the vector sums and
   a transformation matrix that expresses the basis in terms of the 
   vector sums */
			base(&comb[i], &sumvecbase, sumveclist, F.A[0], dim);
			if (flags.PRINT == 1)
/* if the -P option is given, print the rank of the lattice generated by the
   vector sums on level i on AUTO.tmp */
			{
				outfile = fopen("AUTO.tmp", "a");
				fprintf(outfile, "comb[%d].rank = %d\n", i, comb[i].rank);
				fclose(outfile);
			}
/* compute the coefficients of the vector sums in terms of the basis */
			coef(&comb[i], sumvecbase, sumveclist, F.A[0], dim);
			for (j = 0; j <= comb[i].list.n; ++j)
				free(sumveclist[j]);
			free(sumveclist);
/* compute the scalar products of the base-vectors */
			scpforms(&comb[i], sumvecbase, F);
			for (j = 0; j < comb[i].rank; ++j)
				free(sumvecbase[j]);
			free(sumvecbase);
		}
	}
/* the Bacher-polynomials for the first BACHDEP base-vectors are computed,
   if no scalar product was given as an argument, the scalar product is set
   to 1/2 the norm of the base-vector (with respect to the first form) */
	if (flags.BACH[0] == 1)
	{
		bach = (bachpol*)xmalloc(flags.BACHDEP * sizeof(bachpol));
		for (i = 0; i < flags.BACHDEP; ++i)
		{
			if (flags.BACH[2] == 0)
/* no scalar product was given as an option */
				flags.BACHSCP = V.v[fp.e[i]][dim] / 2;
/* compute the Bacher-polynomial */
			bacher(&bach[i], fp.e[i], flags.BACHSCP, V, F.v[0]);
			if (flags.PRINT == 1)
/* if the -P option is given, print the Bacher-polynomial on AUTO.tmp */
			{
				outfile = fopen("AUTO.tmp", "a");
				fprintf(outfile, "Bacher-polynomial of base vector fp.e[%d]:\n", i);
				fputbach(outfile, bach[i]);
				fclose(outfile);
			}
		}
	}
/* set up the group: the generators in G.g[i] are matrices that fix the
   first i base-vectors but do not fix fp.e[i] */
	G.g = (int****)xmalloc(dim * sizeof(int***));
/* G.ng is the number of generators in G.g[i] */
	G.ng = (int*)xmalloc(dim * sizeof(int));
/* the first G.nsg[i] generators in G.g[i] are obtained as stabilizer elements
   and are not necessary to generate the group */
	G.nsg = (int*)xmalloc(dim * sizeof(int));
/* G.ord[i] is the orbit length of fp.e[i] under the generators in 
   G.g[i] ... G.g[dim-1] */
	G.ord = (int*)xmalloc(dim * sizeof(int));
	for (i = 0; i < dim; ++i)
	{
		G.g[i] = (int***)xmalloc(1 * sizeof(int**));
		G.ng[i] = 0;
		G.nsg[i] = 0;
	}
	G.dim = dim;
	G.g[0][0] = (int**)xmalloc(dim * sizeof(int*));
/* -Id is always an automorphism */
	for (i = 0; i < dim; ++i)
	{
		G.g[0][0][i] = (int*)xmalloc(dim * sizeof(int));
		for (j = 0; j < dim; ++j)
			G.g[0][0][i][j] = 0;
		G.g[0][0][i][i] = -1;
		G.ng[0] = 1;
	}
/* fill in the generators which were given as input */
	for (i = 0; i < nH; ++i)
	{
		for (j = 0; j < dim  &&  operate(fp.e[j], H[i], V) == fp.e[j]; ++j);
		if (j < dim)
/* the generator is not Id and fixes fp.e[0],...,fp.e[j-1] but not fp.e[j] */
		{
			++G.ng[j];
			G.g[j] = (int***)xrealloc(G.g[j], G.ng[j] * sizeof(int**));
			G.g[j][G.ng[j]-1] = (int**)xmalloc(dim * sizeof(int*));
			for (k = 0; k < dim; ++k)
			{
				G.g[j][G.ng[j]-1][k] = (int*)xmalloc(dim * sizeof(int));
				for (l = 0; l < dim; ++l)
					G.g[j][G.ng[j]-1][k][l] = H[i][k][l];
			}
		}
		for (k = 0; k < dim; ++k)
			free(H[i][k]);
		free(H[i]);
	}
	if (nH > 0)
		free(H);
	nH = 0;
	for (i = 0; i < dim; ++i)
		nH += G.ng[i];
	H = (int***)xmalloc(nH * sizeof(int**));
/* calculate the orbit lengths under the automorphisms known so far */
	for (i = 0; i < dim; ++i)
	{
		if (G.ng[i] > 0)
		{
			nH = 0;
			for (j = i; j < dim; ++j)
			{
				for (k = 0; k < G.ng[j]; ++k)
				{
					H[nH] = G.g[j][k];
					++nH;
				}
			}
			G.ord[i] = orbitlen(fp.e[i], fp.diag[i], H, nH, V);
		}
		else
			G.ord[i] = 1;
	}
	free(H);
/* compute the full automorphism group */
	autom(&G, V, F, fp, comb, bach, flags);
/* print the generators of the group, which are necessary for generating */
	B = putgens(G, flags);
	n = 0;
	for (i = flags.STAB; i < dim; ++i)
	{
		if (G.ord[i] > 1)
			++n;
	}
/* print a base to AUTO.tmp if the print-flag is set */
	if (flags.PRINT == 1)
	{
		outfile = fopen("AUTO.tmp", "a");
		fprintf(outfile, "%dx%d\n", n, dim);
		for (i = flags.STAB; i < dim; ++i)
		{
			if (G.ord[i] > 1)
			{
				for (j = 0; j < dim; ++j)
					fprintf(outfile, " %d", V.v[fp.e[i]][j]);
				fprintf(outfile, "\n");
			}
		}
		fclose(outfile);
	}
/* print the order of the group */
	putord(G, flags, B);


	if (flags.BACH[0] == 1)
	{
		for (i = 0; i < flags.BACHDEP; ++i)
                   free(bach[i].coef);
                free(bach);
        }
	if (flags.DEPTH > 0)
	{
            for(i=0;i<dim;i++)
            {
              for(j=0;j<=comb[i].list.n;j++)
                free(comb[i].coef[j]);
              free(comb[i].coef);
              for(j=0;j<=comb[i].list.n;j++)
                free(comb[i].list.v[j]);
              free(comb[i].list.v);
              for(j=0;j<comb[i].rank;j++)
                free(comb[i].trans[j]);
              free(comb[i].trans);
              for(j=0;j<F.n;j++)
              {
                 for(k=0;k<comb[i].rank;k++)
                    free(comb[i].F[j][k]);
                 free(comb[i].F[j]);
              }
              free(comb[i].F);
            }
            free(comb);
        }
        for(i=0;i<F.n;i++)
        {
          for(j=1;j<=V.n;j++)
            free(F.v[i][j]);
          free(F.v[i]);
          for(j=0;j<dim;j++)
            free(F.A[i][j]);
          free(F.A[i]);
        }
        free(F.A); free(F.v);
        for(i=0;i<=V.n;i++)
          free(V.v[i]);
        free(V.v);
        free(fp.diag); free(fp.per); free(fp.e);
        for(i=0;i<dim;i++)
        {
           for(j=0;j<G.ng[i];j++)
           {
             for(k=0;k<dim;k++)
               free(G.g[i][j][k]);
             free(G.g[i][j]);
           }
           free(G.g[i]);
        }
        free(G.g); free(G.ord); free(G.ng); free(G.nsg);

        return(B);
}



bravais_TYP *
perfect_normalizer (matrix_TYP *Fo, matrix_TYP *SV, matrix_TYP **Erz, int Erzanz, int *options, matrix_TYP **P, int Panz, int Pdim, matrix_TYP **Pbase)
{
	FILE		*outfile;
	bachpol		*bach;
	flagstruct	flags;
	scpcomb		*comb;
	group		G;
	invar		F;
	veclist		V, norm;
	fpstruct	fp;
	int		dim, max, fail, *vec;
	int		***H, nH, ngen, **sumveclist, **sumvecbase;
	int		i, j, k, l, n, *Vvi, *Vvj, **Fvi, *Fvij, **FAi;
        bravais_TYP *B;

        normal_option = TRUE;
        perp_no = Panz;
        perp = (int ***)xmalloc(perp_no *sizeof(int **));
        for(i=0;i<perp_no;i++)
          perp[i] = P[i]->array.SZ;
        perpdim = Pdim;
        perpbase = (int ***)xmalloc(perpdim *sizeof(int **));
        for(i=0;i<perpdim;i++)
          perpbase[i] = Pbase[i]->array.SZ;
        perpprod = (int ***)xmalloc(perpdim *sizeof(int **));
        for(i=0;i<perpdim;i++)
        {
          perpprod[i] = (int **)xmalloc(Fo->cols *sizeof(int *));
          for(j=0;j<Fo->cols;j++)
          {
            perpprod[i][j] = (int *)xmalloc((j+1) *sizeof(int));
          }
        }
        perpvec = (int *)xmalloc(perpdim *sizeof(int));
        for(i=0;i<perpdim;i++)
          perpprod[i] = init_mat(Fo->cols, Fo->cols, "");
/* get the flags from the command line */
	getflags(&flags, options);
   if(Erzanz > 0)
     flags.GEN = 1;
   else
     flags.GEN = 0;
/* F.n is the number of invariant forms */
	F.n = 1;
	F.A = (int***)xmalloc(1 * sizeof(int**));
/* read the invariant forms */
        F.dim = Fo->cols;
	dim = F.dim;
         if(Fo->cols != F.dim || Fo->rows != F.dim)
         {
            printf("error form 'Fo' in autgrp must be a square matrix\n");
            exit(2);
         }
	  F.A[0] = (int**)xmalloc(dim * sizeof(int*));
          for(j=0;j<F.dim;j++)
          {
	    F.A[0][j] = (int*)xmalloc(dim * sizeof(int));
            for(k=0;k<F.dim;k++)
              F.A[0][j][k] = Fo->array.SZ[j][k];
          }
	if (flags.GEN == 1)
/* some automorphisms are already known */
	{
		ngen = Erzanz;
		H = (int***)xmalloc(ngen * sizeof(int**));
		nH = 0;
		fail = 0;
                for(i=0;i<ngen;i++)
                {
                  if(Erz[i]->cols != dim || Erz[i]->rows != dim)
                  {
                    printf("Elements 'Erz' have wrong dimension\n");
                    exit(3);
                  }
                }
                i = 0;
		while (nH+fail < ngen)
		{
		    H[nH] = (int**)xmalloc(dim * sizeof(int*));
                    for(j=0;j<dim;j++)
                    {
		      H[nH][j] = (int*)xmalloc(dim * sizeof(int));
                      for(k=0;k<dim;k++)
                        H[nH][j][k] = Erz[i]->array.SZ[j][k];
                    }
/* check whether the matrix is really an automorphism, i.e. fixes the forms */
			if (checkgen(H[nH], F) == 0)
/* the matrix is not an automorphism */
			{
				++fail;
				for (j = 0; j < dim; ++j)
					free(H[nH][j]);
                                free(H[nH]);
			}
			else
/* the matrix fixes all forms in F */
				++nH;
                    i++;
		}
	}
	else
		nH = 0;



	V.dim = dim;
	V.len = dim + F.n;
/* get the short vectors */
        V.n = SV->rows;
        V.v = (int**)xmalloc((V.n+1) * sizeof(int*));
        V.v[0] = (int*)xmalloc(V.len * sizeof(int));
        for(i=1;i<= V.n;i++)
        {
               V.v[i] = SV->array.SZ[i-1];
        }
/* the first nonzero entry of each vector is made positive and the maximal
   entry of the vectors is determined */
	max = 1;
	for (i = 1; i <= V.n; ++i)
	{
		Vvi = V.v[i];
		for (j = 0; j < dim  &&  Vvi[j] == 0; ++j);
		if (j < dim  &&  Vvi[j] < 0)
		{
			for (k = j; k < dim; ++k)
				Vvi[k] *= -1;
		}
		for (k = j; k < dim; ++k)
		{
			if (abs(Vvi[k]) > max)
				max = abs(Vvi[k]);
		}
	}
	if (max > MAXENTRY)
/* there is an entry, which is bigger than MAXENTRY, which might cause overflow,
   since the program doesn't use long integer arithmetic */
	{
		fprintf(stderr, "Error: Found entry %d in short vectors.\nTo avoid overflow, the entries should not exceed %d.\n", max, MAXENTRY);
		exit (3);
	}
/* V.prime is a prime p, such that every entry x in the short vectors remains
   unchanged under reduction mod p in symmetric form, i.e. -p/2 < x <= p/2 */
	for (V.prime = 2*max + 1; isprime(V.prime) == 0; ++V.prime);
/* sort the vectors and delete doublets */
	sortvecs(&V);
/* the norm-vector (i.e. the vector of the norms with respect to the different
   invariant forms) of each vector must be equal to the norm-vector of one
   of the vectors from the standard-base */
	norm.v = (int**)xmalloc((dim+1) * sizeof(int*));
	norm.n = dim;
	norm.dim = F.n;
	norm.len = F.n;
	for (i = 1; i <= norm.n; ++i)
	{
/* norm.v[i] is the norm combination of the (i-1)-th base-vector */
		norm.v[i] = (int*)xmalloc(norm.dim * sizeof(int));
		for (k = 0; k < norm.dim; ++k)
			norm.v[i][k] = F.A[k][i-1][i-1];
	}
	sortvecs(&norm);
/* delete those vectors, which can not be the image of any of the standard-base
   vectors */
	checkvecs(&V, F, norm);
	for (i = 1; i <= norm.n; ++i)
		free(norm.v[i]);
	free(norm.v);
/* print the invariant forms, if output is in ASCII-format */
/*********************************************
		printf("#%d\n", F.n);
		for (i = 0; i < F.n; ++i)
			putform(F.A[i], dim);
**********************************************/
/* F.v[i][j] is the transposed of the product of the Gram-matrix F.A[i] with the
   transposed of the vector v[j], hence the scalar product of v[j] and v[k] with
   respect to the Gram-matrix F.A[i] can be computed as standard scalar product
   of v[j] and F.v[i][k] */
	F.v = (int***)xmalloc(F.n * sizeof(int**));
/* the product of the maximal entry in the short vectors with the maximal entry
   in F.v[i] should not exceed MAXNORM to avoid overflow */
	max = MAXNORM / max;
	for (i = 0; i < F.n; ++i)
	{
		FAi = F.A[i];
		F.v[i] = (int**)xmalloc((V.n+1) * sizeof(int*));
		Fvi = F.v[i];
		for (j = 1; j <= V.n; ++j)
		{
			Vvj = V.v[j];
			Fvi[j] = (int*)xmalloc(dim * sizeof(int));
			Fvij = Fvi[j];
			for (k = 0; k < dim; ++k)
			{
				Fvij[k] = sscp(FAi[k], Vvj, dim);
				if (abs(Fvij[k]) > max)
/* some entry in F.v[i] is too large */
				{
					fprintf(stderr, "Error: Found entry %d in F.v.\nTo avoid overflow, the entries should not exceed %d.\n", Fvij[k], max);
					exit (3);
				}
			}
		}
	}
/* fp.per is the order in which the images of the base-vectors are set */
	fp.per = (int*)xmalloc(dim * sizeof(int));
/* fp.e[i] is the index in V.v of the i-th vector of the standard-base in the 
   new order */
	fp.e = (int*)xmalloc(dim * sizeof(int));
/* fp.diag is the diagonal of the fingerprint in the new order, fp.diag[i] is
   an upper bound for the length of the orbit of the i-th base-vector under
   the stabilizer of the preceding ones */
	fp.diag = (int*)xmalloc(dim * sizeof(int));
/* compute the fingerprint */
	fingerprint(&fp, F, V);
	if (flags.PRINT == 1)
/* if the -P option is given, print the diagonal fp.diag and the number of short
   vectors on AUTO.tmp */
	{
		outfile = fopen("AUTO.tmp", "a");
		fprintf(outfile, "%d short vectors\n", 2*V.n);
		fprintf(outfile, "fingerprint diagonal:\n");
		for (i = 0; i < dim; ++i)
			fprintf(outfile, " %2d", fp.diag[i]);
		fprintf(outfile, "\n");
		fclose(outfile);
	}
/* get the standard basis in the new order */
	vec = (int*)xmalloc(dim * sizeof(int));
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
			vec[j] = 0;
		vec[fp.per[i]] = 1;
		fp.e[i] = numberof(vec, V);
		if (abs(fp.e[i]) > V.n)
/* the standard-base must occur in the set of short vectors, since Id is
   certainly an automorphism */
		{
			fprintf(stderr, "Error: standard basis vector nr. %d not found\n", i);
			exit (3);
		}
	}
	free(vec);
/* if the -D option is given, the scalar product combinations and the 
   corresponding vector sums are computed for the standard-basis */
	if (flags.DEPTH > 0)
	{
		comb = (scpcomb*)xmalloc(dim * sizeof(scpcomb));
		for (i = 0; i < dim; ++i)
		{
/* compute the list of scalar product combinations and the corresponding
   vector sums */
			scpvecs(&comb[i].list, &sumveclist, i, fp.e, flags.DEPTH, V, F);
/* compute a basis for the lattice that is generated by the vector sums and
   a transformation matrix that expresses the basis in terms of the 
   vector sums */
			base(&comb[i], &sumvecbase, sumveclist, F.A[0], dim);
			if (flags.PRINT == 1)
/* if the -P option is given, print the rank of the lattice generated by the
   vector sums on level i on AUTO.tmp */
			{
				outfile = fopen("AUTO.tmp", "a");
				fprintf(outfile, "comb[%d].rank = %d\n", i, comb[i].rank);
				fclose(outfile);
			}
/* compute the coefficients of the vector sums in terms of the basis */
			coef(&comb[i], sumvecbase, sumveclist, F.A[0], dim);
			for (j = 0; j <= comb[i].list.n; ++j)
				free(sumveclist[j]);
			free(sumveclist);
/* compute the scalar products of the base-vectors */
			scpforms(&comb[i], sumvecbase, F);
			for (j = 0; j < comb[i].rank; ++j)
				free(sumvecbase[j]);
			free(sumvecbase);
		}
	}
/* the Bacher-polynomials for the first BACHDEP base-vectors are computed,
   if no scalar product was given as an argument, the scalar product is set
   to 1/2 the norm of the base-vector (with respect to the first form) */
	if (flags.BACH[0] == 1)
	{
		bach = (bachpol*)xmalloc(flags.BACHDEP * sizeof(bachpol));
		for (i = 0; i < flags.BACHDEP; ++i)
		{
			if (flags.BACH[2] == 0)
/* no scalar product was given as an option */
				flags.BACHSCP = V.v[fp.e[i]][dim] / 2;
/* compute the Bacher-polynomial */
			bacher(&bach[i], fp.e[i], flags.BACHSCP, V, F.v[0]);
			if (flags.PRINT == 1)
/* if the -P option is given, print the Bacher-polynomial on AUTO.tmp */
			{
				outfile = fopen("AUTO.tmp", "a");
				fprintf(outfile, "Bacher-polynomial of base vector fp.e[%d]:\n", i);
				fputbach(outfile, bach[i]);
				fclose(outfile);
			}
		}
	}
/* set up the group: the generators in G.g[i] are matrices that fix the
   first i base-vectors but do not fix fp.e[i] */
	G.g = (int****)xmalloc(dim * sizeof(int***));
/* G.ng is the number of generators in G.g[i] */
	G.ng = (int*)xmalloc(dim * sizeof(int));
/* the first G.nsg[i] generators in G.g[i] are obtained as stabilizer elements
   and are not necessary to generate the group */
	G.nsg = (int*)xmalloc(dim * sizeof(int));
/* G.ord[i] is the orbit length of fp.e[i] under the generators in 
   G.g[i] ... G.g[dim-1] */
	G.ord = (int*)xmalloc(dim * sizeof(int));
	for (i = 0; i < dim; ++i)
	{
		G.g[i] = (int***)xmalloc(1 * sizeof(int**));
		G.ng[i] = 0;
		G.nsg[i] = 0;
	}
	G.dim = dim;
	G.g[0][0] = (int**)xmalloc(dim * sizeof(int*));
/* -Id is always an automorphism */
	for (i = 0; i < dim; ++i)
	{
		G.g[0][0][i] = (int*)xmalloc(dim * sizeof(int));
		for (j = 0; j < dim; ++j)
			G.g[0][0][i][j] = 0;
		G.g[0][0][i][i] = -1;
		G.ng[0] = 1;
	}
/* fill in the generators which were given as input */
	for (i = 0; i < nH; ++i)
	{
		for (j = 0; j < dim  &&  operate(fp.e[j], H[i], V) == fp.e[j]; ++j);
		if (j < dim)
/* the generator is not Id and fixes fp.e[0],...,fp.e[j-1] but not fp.e[j] */
		{
			++G.ng[j];
			G.g[j] = (int***)xrealloc(G.g[j], G.ng[j] * sizeof(int**));
			G.g[j][G.ng[j]-1] = (int**)xmalloc(dim * sizeof(int*));
			for (k = 0; k < dim; ++k)
			{
				G.g[j][G.ng[j]-1][k] = (int*)xmalloc(dim * sizeof(int));
				for (l = 0; l < dim; ++l)
					G.g[j][G.ng[j]-1][k][l] = H[i][k][l];
			}
		}
		for (k = 0; k < dim; ++k)
			free(H[i][k]);
		free(H[i]);
	}
	if (nH > 0)
		free(H);
	nH = 0;
	for (i = 0; i < dim; ++i)
		nH += G.ng[i];
	H = (int***)xmalloc(nH * sizeof(int**));
/* calculate the orbit lengths under the automorphisms known so far */
	for (i = 0; i < dim; ++i)
	{
		if (G.ng[i] > 0)
		{
			nH = 0;
			for (j = i; j < dim; ++j)
			{
				for (k = 0; k < G.ng[j]; ++k)
				{
					H[nH] = G.g[j][k];
					++nH;
				}
			}
			G.ord[i] = orbitlen(fp.e[i], fp.diag[i], H, nH, V);
		}
		else
			G.ord[i] = 1;
	}
	free(H);
/* compute the full automorphism group */
	autom(&G, V, F, fp, comb, bach, flags);
/* print the generators of the group, which are necessary for generating */
	B = putgens(G, flags);
	n = 0;
	for (i = flags.STAB; i < dim; ++i)
	{
		if (G.ord[i] > 1)
			++n;
	}
/* print a base to AUTO.tmp if the print-flag is set */
	if (flags.PRINT == 1)
	{
		outfile = fopen("AUTO.tmp", "a");
		fprintf(outfile, "%dx%d\n", n, dim);
		for (i = flags.STAB; i < dim; ++i)
		{
			if (G.ord[i] > 1)
			{
				for (j = 0; j < dim; ++j)
					fprintf(outfile, " %d", V.v[fp.e[i]][j]);
				fprintf(outfile, "\n");
			}
		}
		fclose(outfile);
	}
/* print the order of the group */
	putord(G, flags, B);


	if (flags.BACH[0] == 1)
	{
		for (i = 0; i < flags.BACHDEP; ++i)
                   free(bach[i].coef);
                free(bach);
        }
	if (flags.DEPTH > 0)
	{
            for(i=0;i<dim;i++)
            {
              for(j=0;j<=comb[i].list.n;j++)
                free(comb[i].coef[j]);
              free(comb[i].coef);
              for(j=0;j<=comb[i].list.n;j++)
                free(comb[i].list.v[j]);
              free(comb[i].list.v);
              for(j=0;j<comb[i].rank;j++)
                free(comb[i].trans[j]);
              free(comb[i].trans);
              for(j=0;j<F.n;j++)
              {
                 for(k=0;k<comb[i].rank;k++)
                    free(comb[i].F[j][k]);
                 free(comb[i].F[j]);
              }
              free(comb[i].F);
            }
            free(comb);
        }
        for(i=0;i<F.n;i++)
        {
          for(j=1;j<=V.n;j++)
            free(F.v[i][j]);
          free(F.v[i]);
          for(j=0;j<dim;j++)
            free(F.A[i][j]);
          free(F.A[i]);
        }
        free(F.A); free(F.v);
        free(V.v);
        free(perp);
        free(perpbase);
        for(i=0;i<perpdim;i++)
        {
          for(j=0;j<Fo->cols;j++)
            free(perpprod[i][j]);
          free(perpprod[i]);
        }
        free(perpprod);
        free(perpvec);
        free(fp.diag); free(fp.per); free(fp.e);
        for(i=0;i<dim;i++)
        {
           for(j=0;j<G.ng[i];j++)
           {
             for(k=0;k<dim;k++)
               free(G.g[i][j][k]);
             free(G.g[i][j]);
           }
           free(G.g[i]);
        }
        free(G.g); free(G.ord); free(G.ng); free(G.nsg);

        return(B);
}
