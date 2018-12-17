/*****	Main program for the isometry program ISOM	*****/
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: isometry.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

#include"typedef.h"
#include "types.h"

static int normal_option;
static int perp_no;
static int ***perp;
static int perpdim;
static int ***perpbase;
static int ***perpprod;
static int *perpvec;

#include "isotools.c"
#include "bachtools.c"
#include "iotools.c"
#include "lattools.c"
#include "mattools.c"
#include "orbtools.c"
#include "preproc.c"
#include "sorttools.c"
#include "perfecttools.c"
#include "matrix.h"
#include "sort.h"

/*
@-------------------------------------------------------------------------
@ matrix_TYP *isometry(F1, F2, Fanz, SV1, SV2, Erz, Erzanz, options)
@ matrix_TYP **F1, **F2, *SV1, *SV2, **Erz;
@ int Fanz, Erzanz, *options;
@
@ The function 'isometry' calculates a matrix X such that
@    X * F1[i] *X^{tr} = F2[i] for 1<= i<= Foanz
@ returned via a pointer to 'matrix_TYP'.
@ If no such matrix exists the functions returns NULL
@
@ The arguments of isometry are:
@ matrix_TYP	**F1:		a set of n times n matrices,
@				the first must be positiv definite
@ matrix_TYP	**F2:		a set of n times n matrices,
@				the first must be positiv definite
@ int	 	Fanz:		the number of the matrices given in 'F1'.
@ matrix_TYP	*SV1:		The rows of the matrix 'SV1' must be the vectors
@				x in Z^n with x * F1[0] * x^{tr} <= m, where
@				m is the maximal diagonal entry of the Matrix
@				F1[0]
@ matrix_TYP	*SV2:		The rows of the matrix 'SV2' must be the vectors
@				x in Z^n with x * F2[0] * x^{tr} <= m, where
@				m is the maximal diagonal entry of the Matrix
@				F.[0]
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

matrix_TYP *isometry(F1, F2, Fanz, SV1, SV2, Erz, Erzanz, options)
matrix_TYP **F1, **F2, *SV1, *SV2, **Erz;
int Fanz, Erzanz, *options;
{
	FILE		*outfile;
	bachpol		*bach;
	flagstruct	flags;
	scpcomb		*comb;
	invar		F, FF;
	veclist		V, norm;
	fpstruct	fp;
	int		dim, max, fail, *vec;
	int		***G, nG, ngen, **sumveclist, **sumvecbase;
	int		i, j, k, n, *Vvi, *Vvj, **Fvi, *Fvij, **FAi, nV1;
	matrix_TYP	*erg, *erg1;

        extern matrix_TYP *mat_inv();

       normal_option = FALSE;
       if(SV1->rows != SV2->rows)
          return(NULL);
/* FF.v[i][j] is the transposed of the product of the Gram-matrix FF.A[i] 
   get the flags from the command line */
	getflags(&flags, options);
        if(Erzanz > 0)
          flags.GEN = 1;
        else flags.GEN = 0;
/* F.n is the number of invariant forms */
	F.n = Fanz;
	if ((F.A = (int***)malloc(F.n * sizeof(int**))) == 0)
		exit (2);
/* read the invariant forms of the first lattice */
        F.dim = F1[0]->cols;
        dim = F.dim;
        if(F1[0]->cols != F1[0]->rows)
        {
           printf("Non-square matrix as input for isometry\n");
           exit(3);
        }
        for(i=1;i<F.n;i++)
        {
          if(F1[i]->rows != F.dim || F1[i]->cols != F.dim)
          {
             printf("Different dimensions of the forms in isometry\n");
             exit(3);
          }
        }
        for(i=0;i<F.n;i++)
        {
          if((F.A[i] = (int **)malloc(dim *sizeof(int *))) == 0)
          {
             printf("realloc of F.A[%d] failed in isometry\n", i);
             exit(2);
          }
          for(j=0;j<dim;j++)
          {
            if((F.A[i][j] = (int *)malloc(dim *sizeof(int))) == 0)
            {
               printf("realloc of F.A[%d][%d] failed in isometry\n", i, j);
               exit(2);
            }
            for(k=0;k<dim;k++)
             F.A[i][j][k] = F1[i]->array.SZ[j][k];
          }
        }
/* get the short vectors of the first lattice */
	V.dim = dim;
	V.len = dim + F.n;
        V.n = SV1->rows;
        if((V.v = (int **)malloc((V.n+1) *sizeof(int *))) == 0)
        {
           printf("realloc of V.v failed in isometry\n");
           exit(2);
        }
        if((V.v[0] = (int *)malloc(V.len *sizeof(int))) == 0)
        {
           printf("realloc of V.v[0] failed in isometry\n");
             exit(2);
        }
        for(j=0;j<V.len;j++)
           V.v[0][j] = 0;
        for(i=1;i<=V.n;i++)
        {
          if((V.v[i] = (int *)malloc(V.len *sizeof(int))) == 0)
          {
             printf("realloc of V.v[%d] failed in isometry\n", i);
             exit(2);
          }
          for(j=0;j<=V.dim;j++)
             V.v[i][j] = SV1->array.SZ[i-1][j];
        }
/* the first nonzero entry of each vector is made positive and the maximal 
   entry in the vectors is determined */
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
		fprintf(stderr, "Error: Found entry %d in short vectors of the first lattice.\nTo avoid overflow, the entries should not exceed %d.\n", max, MAXENTRY);
		exit (3);
	}
/* sort the vectors and delete doublets */
	sortvecs(&V);
/* the norm-vector (i.e. the vector of the norms with respect to the different
   invariant forms) of each vector must be equal to the norm-vector of one
   of the vectors from the standard-base */
	if ((norm.v = (int**)malloc((dim+1) * sizeof(int*))) == 0)
		exit (2);
	norm.n = dim;
	norm.dim = F.n;
	norm.len = F.n;
	for (i = 1; i <= norm.n; ++i)
	{
/* norm.v[i] is the norm combination of the (i-1)-th base-vector */
		if ((norm.v[i] = (int*)malloc(norm.dim * sizeof(int))) == 0)
			exit (2);
		for (k = 0; k < norm.dim; ++k)
			norm.v[i][k] = F.A[k][i-1][i-1];
	}
	sortvecs(&norm);
/* delete those vectors, which can not be the image of any vector of the
   standard-base */
	checkvecs(&V, F, norm);
/* F.v[i][j] is the transposed of the product of the Gram-matrix F.A[i] with the
   transposed of the vector v[j], hence the scalar product of v[j] and v[k] with
   respect to the Gram-matrix F.A[i] can be computed as standard scalar product
   of v[j] and F.v[i][k] */
	if ((F.v = (int***)malloc(F.n * sizeof(int**))) == 0)
		exit (2);
/* the product of the maximal entry in the short vectors with the maximal entry
   in F.v[i] should not exceed MAXNORM to avoid overflow */
	max = MAXNORM / max;
	for (i = 0; i < F.n; ++i)
	{
		FAi = F.A[i];
		if ((F.v[i] = (int**)malloc((V.n+1) * sizeof(int*))) == 0)
			exit (2);
		Fvi = F.v[i];
		for (j = 1; j <= V.n; ++j)
		{
			Vvj = V.v[j];
			if ((Fvi[j] = (int*)malloc(dim * sizeof(int))) == 0)
				exit (2);
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
/* fp.per is the order in which the images of the base vectors are chosen */
	if ((fp.per = (int*)malloc(dim * sizeof(int))) == 0)
		exit (2);
/* fp.e[i] is the index in V.v of the i-th vector of the standard-base in the 
   new order */
	if ((fp.e = (int*)malloc(dim * sizeof(int))) == 0)
		exit (2);
/* fp.diag is the diagonal of the fingerprint in the new order */
	if ((fp.diag = (int*)malloc(dim * sizeof(int))) == 0)
		exit (2);
/* compute the fingerprint */
	fingerprint(&fp, F, V);
	if (flags.PRINT == 1)
/* if the -P option is given, print the diagonal fp.diag and the number of
   short vectors on ISOM.tmp */
	{
		outfile = fopen("ISOM.tmp", "a");
		fprintf(outfile, "%d short vectors\n", 2*V.n);
		fprintf(outfile, "fingerprint diagonal:\n");
		for (i = 0; i < dim; ++i)
			fprintf(outfile, " %2d", fp.diag[i]);
		fprintf(outfile, "\n");
		fclose(outfile);
	}
/* get the standard basis in the new order */
	if ((vec = (int*)malloc(dim * sizeof(int))) == 0)
		exit (2);
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
			vec[j] = 0;
		vec[fp.per[i]] = 1;
		fp.e[i] = numberof(vec, V);
		if (abs(fp.e[i]) > V.n)
		{
			fprintf(stderr, "Error: standard basis vector nr. %d not found\n", i);
			exit (3);
		}
	}
	free(vec);
/* if the -D option is given, the scalar product combinations and the 
   corresponding vector sums are computed for the basis of the first lattice */
	if (flags.DEPTH > 0)
	{
		if ((comb = (scpcomb*)malloc(dim * sizeof(scpcomb))) == 0)
			exit (2);
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
   vector sums on level i on ISOM.tmp */
			{
				outfile = fopen("ISOM.tmp", "a");
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
/* the Bacher-polynomials for the first BACHDEP base vectors are computed,
   if no scalar product was given as an argument, the scalar product is set
   to 1/2 the norm of the base-vector (with respect to the first form) */
	if (flags.BACH[0] == 1)
	{
		if ((bach = (bachpol*)malloc(flags.BACHDEP * sizeof(bachpol))) == 0)
			exit (2);
		for (i = 0; i < flags.BACHDEP; ++i)
		{
			if (flags.BACH[2] == 0)
/* no scalar product was given as an option */
				flags.BACHSCP = V.v[fp.e[i]][dim] / 2;
/* compute the Bacher-polynomial */
			bacher(&bach[i], fp.e[i], flags.BACHSCP, V, F.v[0]);
			if (flags.PRINT == 1)
/* if the -P option is given, print the Bacher-polynomial on ISOM.tmp */
			{
				outfile = fopen("ISOM.tmp", "a");
				fprintf(outfile, "Bacher-polynomial of base vector fp.e[%d]:\n", i);
				fputbach(outfile, bach[i]);
				fclose(outfile);
			}
		}
	}
/* the vectors of the first lattice are no longer required, only their number */
	nV1 = V.n;
        for(i=0;i<=V.n;i++)
          free(V.v[i]);
        free(V.v);
	for (i = 0; i < F.n; ++i)
	{
		for (j = 1; j <= V.n; ++j)
			free(F.v[i][j]);
		free(F.v[i]);
	}
	free(F.v);
/* FF are the invariant forms of the second lattice */
	if ((FF.A = (int***)malloc(F.n * sizeof(int**))) == 0)
		exit (2);
	FF.n = F.n;
	FF.dim = dim;
/* read the invariant forms of the second lattice */
        for(i=1;i<F.n;i++)
        {
          if(F2[i]->rows != F.dim || F2[i]->cols != F.dim)
          {
             printf("Different dimensions of the forms in isometry\n");
             exit(2);
          }
        }
        for(i=0;i<F.n;i++)
        {
          if((FF.A[i] = (int **)malloc(dim *sizeof(int *))) == 0)
          {
             printf("realloc of FF.A[%d] failed in isometry\n", i);
             exit(2);
          }
          for(j=0;j<dim;j++)
          {
            if((FF.A[i][j] = (int *)malloc(dim *sizeof(int))) == 0)
            {
               printf("realloc of FF.A[%d][%d] failed in isometry\n", i, j);
               exit(2);
            }
            for(k=0;k<dim;k++)
             FF.A[i][j][k] = F2[i]->array.SZ[j][k];
          }
        }
/* get the short vectors of the second lattice */
	V.dim = dim;
	V.len = dim + F.n;
        V.n = SV2->rows;
        if((V.v = (int **)malloc((V.n+1) *sizeof(int *))) == 0)
        {
           printf("realloc of V.v failed in isometry\n");
           exit(2);
        }
        if((V.v[0] = (int *)malloc(V.len *sizeof(int))) == 0)
        {
           printf("realloc of V.v[0] failed in isometry\n");
             exit(2);
        }
        for(j=0;j<V.len;j++)
           V.v[0][j] = 0;
        for(i=1;i<=V.n;i++)
        {
          if((V.v[i] = (int *)malloc(V.len *sizeof(int))) == 0)
          {
             printf("realloc of V.v[%d] failed in isometry\n", i);
             exit(2);
          }
          for(j=0;j<=V.dim;j++)
           V.v[i][j] = SV2->array.SZ[i-1][j];
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
		fprintf(stderr, "Error: Found entry %d in short vectors of the second lattice.\nTo avoid overflow, the entries should not exceed %d.\n", max, MAXENTRY);
		exit (3);
	}
/* V.prime is a prime p, such that the entries of the short vectors remain 
   unchanged under reduction mod p in symmetric form, i.e. -p/2 < x <= p/2 */
	for (V.prime = 2*max + 1; isprime(V.prime) == 0; ++V.prime);
        j = F.A[0][0][0];
        for(k=1;k<dim;k++)
        {
           if(F.A[0][k][k] > j)
             j = F.A[0][k][k];
        }
        for(k=0;k<dim;k++)
        {
           if(FF.A[0][k][k] > j)
             V.prime = DEF_PRIME;
        }
/* sort the vectors and delete doublets */
	sortvecs(&V);
/* the norm-vector (i.e. the vector of the norms with respect to the different
   invariant forms) of each vector must be equal to the norm-vector of one
   of the vectors from the standard-base of the first lattice */
	checkvecs(&V, FF, norm);
	for (i = 1; i <= norm.n; ++i)
		free(norm.v[i]);
	free(norm.v);
/* FF.v[i][j] is the transposed of the product of the Gram-matrix FF.A[i] 
   with the transposed of the vector V.v[j], which now is a short vector of
   the second lattice */
	if ((FF.v = (int***)malloc(FF.n * sizeof(int**))) == 0)
		exit (2);
/* the product of the maximal entry in the short vectors with the maximal entry
   in FF.v[i] should not exceed MAXNORM to avoid overflow */
	max = MAXNORM / max;
	for (i = 0; i < FF.n; ++i)
	{
		FAi = FF.A[i];
		if ((FF.v[i] = (int**)malloc((V.n+1) * sizeof(int*))) == 0)
			exit (2);
		Fvi = FF.v[i];
		for (j = 1; j <= V.n; ++j)
		{
			Vvj = V.v[j];
			if ((Fvi[j] = (int*)malloc(dim * sizeof(int))) == 0)
				exit (2);
			Fvij = Fvi[j];
			for (k = 0; k < dim; ++k)
			{
				Fvij[k] = sscp(FAi[k], Vvj, dim);
				if (abs(Fvij[k]) > max)
/* some entry in FF.v[i] is too large */
				{
					fprintf(stderr, "Error: Found entry %d in FF.v.\nTo avoid overflow, the entries should not exceed %d.\n", Fvij[k], max);
					exit (3);
				}
			}
		}
	}
/* G are the automorphisms of the second lattice */
	if ((G = (int***)malloc(1 * sizeof(int**))) == 0)
		exit (2);
/* nG is the number of automorphisms of the second lattice */
	nG = 0;
	if (flags.GEN == 1)
/* get the automorphisms, which are already known */
	{
		ngen = Erzanz;
		if ((G = (int***)realloc(G, ngen * sizeof(int**))) == 0)
			exit (2);
		fail = 0;
                i=0;
		while (nG+fail < ngen)
		{
			if (Erz[i]->cols != dim || Erz[i]->rows != dim)
			{
				fprintf(stderr, "Error: dimension %d should be %d\n", n, dim);
				exit (3);
			}
                        if((G[nG] = (int **)malloc(dim *sizeof(int *))) == 0)
                        {
                          printf("realloc of G[%d] in isometry failed\n", nG); 
                          exit(2);
                        }
                        for(j=0;j<dim;j++)
                        {
                          if((G[nG][j] = (int *)malloc(dim *sizeof(int))) == 0)
                          {
                            printf("realloc of G[%d][%d] in isometry failed\n", nG, j); 
                            exit(2);
                          }
                          for(k=0;k<dim;k++)
                              G[nG][j][k] = Erz[i]->array.SZ[k][j];
                        }
/* check whether the matrix is really an automorphism, i.e. fixes the forms */
			if (checkgen(G[nG], FF) == 0)
/* the matrix is not an automorphism */
			{
				++fail;
				for (j = 0; j < dim; ++j)
					free(G[nG][j]);
                                free(G[nG]);
			}
			else
/* the matrix fixes the forms in FF */
				++nG;
                        i++;
		}
	}
	if (nG == 0)
/* if no automorphisms are known at least take -Id */
	{
		if ((G[0] = (int**)malloc(dim * sizeof(int*))) == 0)
			exit (2);
		for (i = 0; i < dim; ++i)
		{
			if ((G[0][i] = (int*)malloc(dim * sizeof(int))) == 0)
				exit (2);
			for (j = 0; j < dim; ++j)
				G[0][i][j] = 0;
			G[0][i][i] = -1;
		}
		nG = 1;
	}
/* now search for an isometry */
	erg = bs_isometry(V, F, FF, fp, G, nG, comb, bach, flags);
/**********************
        for(i=0;i<nG;i++)
        {
           for(j=0;j<dim;j++)
             free(G[i][j]);
           free(G[i]);
        }
*************************/
        free(G);
        for(i=0;i<F.n;i++)
        {
          for(j=0;j<dim;j++)
            free(F.A[i][j]);
          free(F.A[i]);
        }
        free(F.A);
        for(i=0;i<FF.n;i++)
        {
          for(j=1;j<=V.n;j++)
            free(FF.v[i][j]);
          free(FF.v[i]);
          for(j=0;j<dim;j++)
            free(FF.A[i][j]);
          free(FF.A[i]);
        }
        free(FF.A); free(FF.v);
        for(i=0;i<=V.n;i++)
          free(V.v[i]);
        free(V.v);
	if (flags.BACH[0] == 1)
	{
            for(i=0;i<flags.BACHDEP;i++)
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
        free(fp.diag); free(fp.per); free(fp.e);
        if(erg == NULL)
          return(NULL);
        erg1 = mat_inv(erg);
        free_mat(erg);
        return(erg1);
}

/*
@-------------------------------------------------------------------------
@ matrix_TYP *perfect_normal_isometry(F1, F2, SV1, SV2, Erz, Erzanz,
@                                     options, P, Panz, Pbase, Pdim)
@ matrix_TYP *F1, *F2, *SV1, *SV2, **Erz, **P, **Pbase;
@ int Erzanz, *options, Panz, Pdim;
@
@ The function 'perfect_normal_isometry' calculates a matrix X with
@                X * F1 * X^{tr} = F2 and
@    X^{-1} * Pbase[i] * X^{-tr} is a matrix in P for 0<= i< Pdim
@
@
@ The arguments are:
@ F1:       a symmetric positive definite n by n-matrix.
@ F2:       a symmetric positive definite n by n-matrix.
@ SV1:      the vectors x in GL_n(Z) (up to sign) with xF1x^{tr} < k,
@           where k is the maximal diagonal entry of 'F1'.
@ SV2:      the vectors x in GL_n(Z) (up to sign) with xF2x^{tr} < k,
@           where k is the maximal diagonal entry of 'F1'.
@ Erz:      matrices of elements of Aut(F2) (if known).
@ Erzanz:   the number of matrices in 'Erz'.
@ options:  the same as for  the function 'isometry'.
@ P:        a set of matrices
@ Panz:     the number of matrices given in 'P'.
@ Pbase:    a set of matrices
@ Pdim:     the number of matrices given in 'Pbase'.
@
@ This function is designed to calculate  an isometry of two
@ G-perfect forms that is in the normalizer of a group G (the
@ generators of G are the matrices 'Erz').
@ Then 'P' is the set of matrices given by the Voronoi-directions of 'F1'
@ and 'Pbase' a maximal linearly independent subset of the 
@ Voronoi-directions of 'F2'.
@
@---------------------------------------------------------------------------
@
*/


matrix_TYP *perfect_normal_isometry(F1, F2, SV1, SV2, Erz, Erzanz, options, P, Panz, Pbase, Pdim)
matrix_TYP *F1, *F2, *SV1, *SV2, **Erz, **P, **Pbase;
int Erzanz, *options, Panz, Pdim;
{
	FILE		*outfile;
	bachpol		*bach;
	flagstruct	flags;
	scpcomb		*comb;
	invar		F, FF;
	veclist		V, norm;
	fpstruct	fp;
	int		dim, max, fail, *vec;
	int		***G, nG, ngen, **sumveclist, **sumvecbase;
	int		i, j, k, n, *Vvi, *Vvj, **Fvi, *Fvij, **FAi, nV1;
	matrix_TYP	*erg, *erg1;

        extern int pointer_lower_triangular_mat_comp();
        extern matrix_TYP *mat_inv();

       if(SV1->rows != SV2->rows)
          return(NULL);
       normal_option = TRUE;
       perpdim = Pdim;
       perp_no = Panz;
/* FF.v[i][j] is the transposed of the product of the Gram-matrix FF.A[i] 
   get the flags from the command line */
	getflags(&flags, options);
        if(Erzanz > 0)
          flags.GEN = 1;
        else flags.GEN = 0;
/* F.n is the number of invariant forms */
	F.n = 1;
	if ((F.A = (int***)malloc(F.n * sizeof(int**))) == 0)
		exit (2);
/* read the invariant forms of the first lattice */
        F.dim = F1->cols;
        dim = F.dim;
        if(F1->cols != F1->rows)
        {
           printf("Non-square matrix as input for isometry\n");
           exit(2);
        }
          if((F.A[0] = (int **)malloc(dim *sizeof(int *))) == 0)
          {
             printf("realloc of F.A[0] failed in isometry\n");
             exit(2);
          }
          for(j=0;j<dim;j++)
          {
            if((F.A[0][j] = (int *)malloc(dim *sizeof(int))) == 0)
            {
               printf("realloc of F.A[0][%d] failed in isometry\n", j);
               exit(2);
            }
            for(k=0;k<dim;k++)
             F.A[0][j][k] = F1->array.SZ[j][k];
          }
/* get the short vectors of the first lattice */
	V.dim = dim;
	V.len = dim + F.n;
        V.n = SV1->rows;
        if((V.v = (int **)malloc((V.n+1) *sizeof(int *))) == 0)
        {
           printf("realloc of V.v failed in isometry\n");
           exit(2);
        }
        if((V.v[0] = (int *)malloc(V.len *sizeof(int))) == 0)
        {
           printf("realloc of V.v[0] failed in isometry\n");
             exit(2);
        }
        for(j=0;j<V.len;j++)
           V.v[0][j] = 0;
        for(i=1;i<=V.n;i++)
        {
          if((V.v[i] = (int *)malloc(V.len *sizeof(int))) == 0)
          {
             printf("realloc of V.v[%d] failed in isometry\n", i);
             exit(2);
          }
          for(j=0;j<=V.dim;j++)
             V.v[i][j] = SV1->array.SZ[i-1][j];
        }
/* the first nonzero entry of each vector is made positive and the maximal 
   entry in the vectors is determined */
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
		fprintf(stderr, "Error: Found entry %d in short vectors of the first lattice.\nTo avoid overflow, the entries should not exceed %d.\n", max, MAXENTRY);
		exit (3);
	}
/* sort the vectors and delete doublets */
	sortvecs(&V);
/* the norm-vector (i.e. the vector of the norms with respect to the different
   invariant forms) of each vector must be equal to the norm-vector of one
   of the vectors from the standard-base */
	if ((norm.v = (int**)malloc((dim+1) * sizeof(int*))) == 0)
		exit (2);
	norm.n = dim;
	norm.dim = F.n;
	norm.len = F.n;
	for (i = 1; i <= norm.n; ++i)
	{
/* norm.v[i] is the norm combination of the (i-1)-th base-vector */
		if ((norm.v[i] = (int*)malloc(norm.dim * sizeof(int))) == 0)
			exit (2);
		for (k = 0; k < norm.dim; ++k)
			norm.v[i][k] = F.A[k][i-1][i-1];
	}
	sortvecs(&norm);
/* delete those vectors, which can not be the image of any vector of the
   standard-base */
	checkvecs(&V, F, norm);
/* F.v[i][j] is the transposed of the product of the Gram-matrix F.A[i] with the
   transposed of the vector v[j], hence the scalar product of v[j] and v[k] with
   respect to the Gram-matrix F.A[i] can be computed as standard scalar product
   of v[j] and F.v[i][k] */
	if ((F.v = (int***)malloc(F.n * sizeof(int**))) == 0)
		exit (2);
/* the product of the maximal entry in the short vectors with the maximal entry
   in F.v[i] should not exceed MAXNORM to avoid overflow */
	max = MAXNORM / max;
	for (i = 0; i < F.n; ++i)
	{
		FAi = F.A[i];
		if ((F.v[i] = (int**)malloc((V.n+1) * sizeof(int*))) == 0)
			exit (2);
		Fvi = F.v[i];
		for (j = 1; j <= V.n; ++j)
		{
			Vvj = V.v[j];
			if ((Fvi[j] = (int*)malloc(dim * sizeof(int))) == 0)
				exit (2);
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
/* fp.per is the order in which the images of the base vectors are chosen */
	if ((fp.per = (int*)malloc(dim * sizeof(int))) == 0)
		exit (2);
/* fp.e[i] is the index in V.v of the i-th vector of the standard-base in the 
   new order */
	if ((fp.e = (int*)malloc(dim * sizeof(int))) == 0)
		exit (2);
/* fp.diag is the diagonal of the fingerprint in the new order */
	if ((fp.diag = (int*)malloc(dim * sizeof(int))) == 0)
		exit (2);
/* compute the fingerprint */
	fingerprint(&fp, F, V);
        mach_perp_matrices(fp, P, Pbase, dim);
        pointer_mat_quicksort(perp, 0, perp_no-1, dim, dim, pointer_lower_triangular_mat_comp);
	if (flags.PRINT == 1)
/* if the -P option is given, print the diagonal fp.diag and the number of
   short vectors on ISOM.tmp */
	{
		outfile = fopen("ISOM.tmp", "a");
		fprintf(outfile, "%d short vectors\n", 2*V.n);
		fprintf(outfile, "fingerprint diagonal:\n");
		for (i = 0; i < dim; ++i)
			fprintf(outfile, " %2d", fp.diag[i]);
		fprintf(outfile, "\n");
		fclose(outfile);
	}
/* get the standard basis in the new order */
	if ((vec = (int*)malloc(dim * sizeof(int))) == 0)
		exit (2);
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
			vec[j] = 0;
		vec[fp.per[i]] = 1;
		fp.e[i] = numberof(vec, V);
		if (abs(fp.e[i]) > V.n)
		{
			fprintf(stderr, "Error: standard basis vector nr. %d not found\n", i);
			exit (3);
		}
	}
	free(vec);
/* if the -D option is given, the scalar product combinations and the 
   corresponding vector sums are computed for the basis of the first lattice */
	if (flags.DEPTH > 0)
	{
		if ((comb = (scpcomb*)malloc(dim * sizeof(scpcomb))) == 0)
			exit (2);
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
   vector sums on level i on ISOM.tmp */
			{
				outfile = fopen("ISOM.tmp", "a");
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
/* the Bacher-polynomials for the first BACHDEP base vectors are computed,
   if no scalar product was given as an argument, the scalar product is set
   to 1/2 the norm of the base-vector (with respect to the first form) */
	if (flags.BACH[0] == 1)
	{
		if ((bach = (bachpol*)malloc(flags.BACHDEP * sizeof(bachpol))) == 0)
			exit (2);
		for (i = 0; i < flags.BACHDEP; ++i)
		{
			if (flags.BACH[2] == 0)
/* no scalar product was given as an option */
				flags.BACHSCP = V.v[fp.e[i]][dim] / 2;
/* compute the Bacher-polynomial */
			bacher(&bach[i], fp.e[i], flags.BACHSCP, V, F.v[0]);
			if (flags.PRINT == 1)
/* if the -P option is given, print the Bacher-polynomial on ISOM.tmp */
			{
				outfile = fopen("ISOM.tmp", "a");
				fprintf(outfile, "Bacher-polynomial of base vector fp.e[%d]:\n", i);
				fputbach(outfile, bach[i]);
				fclose(outfile);
			}
		}
	}
/* the vectors of the first lattice are no longer required, only their number */
	nV1 = V.n;
        for(i=0;i<=V.n;i++)
          free(V.v[i]);
        free(V.v);
	for (i = 0; i < F.n; ++i)
	{
		for (j = 1; j <= V.n; ++j)
			free(F.v[i][j]);
		free(F.v[i]);
	}
	free(F.v);
/* FF are the invariant forms of the second lattice */
	if ((FF.A = (int***)malloc(F.n * sizeof(int**))) == 0)
		exit (2);
	FF.n = F.n;
	FF.dim = dim;
/* read the invariant forms of the second lattice */
          if((FF.A[0] = (int **)malloc(dim *sizeof(int *))) == 0)
          {
             printf("realloc of FF.A[0] failed in isometry\n");
             exit(2);
          }
          for(j=0;j<dim;j++)
          {
            if((FF.A[0][j] = (int *)malloc(dim *sizeof(int))) == 0)
            {
               printf("realloc of FF.A[0][%d] failed in isometry\n", j);
               exit(2);
            }
            for(k=0;k<dim;k++)
             FF.A[0][j][k] = F2->array.SZ[j][k];
          }
/* get the short vectors of the second lattice */
	V.dim = dim;
	V.len = dim + F.n;
        V.n = SV2->rows;
        if((V.v = (int **)malloc((V.n+1) *sizeof(int *))) == 0)
        {
           printf("realloc of V.v failed in isometry\n");
           exit(2);
        }
        if((V.v[0] = (int *)malloc(V.len *sizeof(int))) == 0)
        {
           printf("realloc of V.v[0] failed in isometry\n");
             exit(2);
        }
        for(j=0;j<V.len;j++)
           V.v[0][j] = 0;
        for(i=1;i<=V.n;i++)
        {
          if((V.v[i] = (int *)malloc(V.len *sizeof(int))) == 0)
          {
             printf("realloc of V.v[%d] failed in isometry\n", i);
             exit(2);
          }
          for(j=0;j<=V.dim;j++)
             V.v[i][j] = SV2->array.SZ[i-1][j];
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
		fprintf(stderr, "Error: Found entry %d in short vectors of the second lattice.\nTo avoid overflow, the entries should not exceed %d.\n", max, MAXENTRY);
		exit (3);
	}
/* V.prime is a prime p, such that the entries of the short vectors remain 
   unchanged under reduction mod p in symmetric form, i.e. -p/2 < x <= p/2 */
	for (V.prime = 2*max + 1; isprime(V.prime) == 0; ++V.prime);
/* sort the vectors and delete doublets */
	sortvecs(&V);
/* the norm-vector (i.e. the vector of the norms with respect to the different
   invariant forms) of each vector must be equal to the norm-vector of one
   of the vectors from the standard-base of the first lattice */
	checkvecs(&V, FF, norm);
	for (i = 1; i <= norm.n; ++i)
		free(norm.v[i]);
	free(norm.v);
/* FF.v[i][j] is the transposed of the product of the Gram-matrix FF.A[i] 
   with the transposed of the vector V.v[j], which now is a short vector of
   the second lattice */
	if ((FF.v = (int***)malloc(FF.n * sizeof(int**))) == 0)
		exit (2);
/* the product of the maximal entry in the short vectors with the maximal entry
   in FF.v[i] should not exceed MAXNORM to avoid overflow */
	max = MAXNORM / max;
	for (i = 0; i < FF.n; ++i)
	{
		FAi = FF.A[i];
		if ((FF.v[i] = (int**)malloc((V.n+1) * sizeof(int*))) == 0)
			exit (2);
		Fvi = FF.v[i];
		for (j = 1; j <= V.n; ++j)
		{
			Vvj = V.v[j];
			if ((Fvi[j] = (int*)malloc(dim * sizeof(int))) == 0)
				exit (2);
			Fvij = Fvi[j];
			for (k = 0; k < dim; ++k)
			{
				Fvij[k] = sscp(FAi[k], Vvj, dim);
				if (abs(Fvij[k]) > max)
/* some entry in FF.v[i] is too large */
				{
					fprintf(stderr, "Error: Found entry %d in FF.v.\nTo avoid overflow, the entries should not exceed %d.\n", Fvij[k], max);
					exit (3);
				}
			}
		}
	}
/* G are the automorphisms of the second lattice */
	if ((G = (int***)malloc(1 * sizeof(int**))) == 0)
		exit (2);
/* nG is the number of automorphisms of the second lattice */
	nG = 0;
	if (flags.GEN == 1)
/* get the automorphisms, which are already known */
	{
		ngen = Erzanz;
		if ((G = (int***)realloc(G, ngen * sizeof(int**))) == 0)
			exit (2);
		fail = 0;
                i=0;
		while (nG+fail < ngen)
		{
			if (Erz[i]->cols != dim || Erz[i]->rows != dim)
			{
				fprintf(stderr, "Error: dimension %d should be %d\n", n, dim);
				exit (3);
			}
                        if((G[nG] = (int **)malloc(dim *sizeof(int *))) == 0)
                        {
                          printf("realloc of G[%d] in isometry failed\n", nG); 
                          exit(2);
                        }
                        for(j=0;j<dim;j++)
                        {
                          if((G[nG][j] = (int *)malloc(dim *sizeof(int))) == 0)
                          {
                            printf("realloc of G[%d][%d] in isometry failed\n", nG, j); 
                            exit(2);
                          }
                          for(k=0;k<dim;k++)
                              G[nG][j][k] = Erz[i]->array.SZ[k][j];
                        }
/* check whether the matrix is really an automorphism, i.e. fixes the forms */
			if (checkgen(G[nG], FF) == 0)
/* the matrix is not an automorphism */
			{
				++fail;
				for (j = 0; j < dim; ++j)
					free(G[nG][j]);
                                free(G[nG]);
			}
			else
/* the matrix fixes the forms in FF */
				++nG;
                        i++;
		}
	}
	if (nG == 0)
/* if no automorphisms are known at least take -Id */
	{
		if ((G[0] = (int**)malloc(dim * sizeof(int*))) == 0)
			exit (2);
		for (i = 0; i < dim; ++i)
		{
			if ((G[0][i] = (int*)malloc(dim * sizeof(int))) == 0)
				exit (2);
			for (j = 0; j < dim; ++j)
				G[0][i][j] = 0;
			G[0][i][i] = -1;
		}
		nG = 1;
	}
/* now search for an isometry */
	erg = bs_isometry(V, F, FF, fp, G, nG, comb, bach, flags);
        for(i=0;i<F.n;i++)
        {
          for(j=0;j<dim;j++)
            free(F.A[i][j]);
          free(F.A[i]);
        }
        free(F.A);
        for(i=0;i<FF.n;i++)
        {
          for(j=1;j<=V.n;j++)
            free(FF.v[i][j]);
          free(FF.v[i]);
          for(j=0;j<dim;j++)
            free(FF.A[i][j]);
          free(FF.A[i]);
        }
        free(FF.A); free(FF.v);
        for(i=0;i<=V.n;i++)
          free(V.v[i]);
        free(V.v);
	if (flags.BACH[0] == 1)
	{
            for(i=0;i<flags.BACHDEP;i++)
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
        free(fp.diag); free(fp.per); free(fp.e);
/****************************
        for(i=0;i<nG;i++)
        {
           for(j=0;j<dim;j++)
             free(G[i][j]);
           free(G[i]);
        }
*********************/
        free(G);
        free_perp_matrices(dim);
        if(erg == NULL)
          return(NULL);
        erg1 = mat_inv(erg);
        free_mat(erg);
        return(erg1);
}
