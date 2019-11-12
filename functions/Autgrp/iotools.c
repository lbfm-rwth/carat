/*****	This file contains some routines for input/output	*****/
#include "typedef.h"
#include "utils.h"
#include "types.h"
#include "matrix.h"


/*****************************************************\
|	gets the options from the command line
\*****************************************************/
void 
getflags (flagstruct *fl, int *options)
{
 if(options == NULL)
 {
    fl->DEPTH = 0; 
    fl->STAB = 0;
    fl->PRINT = 0;
    fl->BACH[0] = fl->BACH[1] = fl->BACH[2] = 0;
    fl->BACHDEP = 0;
    fl->BACHSCP = 0;
    return;
 }

/* depth for the scalar product combinations */
        if(options[0] >= 0)
          fl->DEPTH = options[0] ;
        else
	  fl->DEPTH = 0;
/* only the point stabilizer of the first STAB basis-vectors will be computed */
        if(options[1] >= 0)
	  fl->STAB = options[1];
        else
	  fl->STAB = 0;
/* flag that every new generator is immediately written on the file AUTO.tmp */
        if(options[2] == 1)
	  fl->PRINT = 1;
        else
	  fl->PRINT = 0;
/* flag that Bacher-polynomials will be used */
        if(options[3]==1 || options[3]==2 || options[3]==3 || options[3]==4)
	  fl->BACH[0] = 1;
        else
	  fl->BACH[0] = 0;
/* flag that the depth for the Bacher-polynomials is given as an argument, 
   default is 1 */
        if(options[3]==2 || options[3]==4)
	  fl->BACH[1] = 1;
        else
	  fl->BACH[1] = 0;
/* flag that the scalar product for the Bacher-polynomials is given as an 
   argument, default is 1/2*norm of the vector */
        if(options[3]==3 || options[3]==4)
	  fl->BACH[2] = 1;
        else
	  fl->BACH[2] = 0;

/* depth for the Bacher-polynomials */
        if(fl->BACH[1] == 1 && options[4] > 0)
          fl->BACHDEP = options[4];
        else
	  fl->BACHDEP = 1;
        if(fl->BACH[0] == 0)
	  fl->BACHDEP = 0;
/* scalar product for the Bacher-polynomials */
        if(fl->BACH[2] == 1)
	   fl->BACHSCP = options[5];
        else
           fl->BACHSCP = 0;
}





/**************************************************************\
|   saves the generators of the group, which are
|   necessary for generating, in a bravais_TYP.
\**************************************************************/
bravais_TYP *
putgens (group G, flagstruct flags)
{
	int	i, j, k, l, dim, ngen, nr;
        bravais_TYP *B;

        B = (bravais_TYP *)malloc(sizeof(bravais_TYP));
        B->gen_no = 0;
        B->form_no = 0;
        B->zentr_no = 0;
        B->normal_no = 0;
        B->cen_no = 0;
        B->dim = G.dim;

	dim = G.dim;
	ngen = 0;
	for (i = flags.STAB; i < dim; ++i)
		ngen += G.ng[i] - G.nsg[i];
        B->gen_no = ngen;
        if(ngen > 0)
        {
           B->gen = (matrix_TYP **)xmalloc(ngen *sizeof(matrix_TYP *));
        }
	nr = 0;
	for (i = flags.STAB; i < dim; ++i)
	{
		for (j = G.nsg[i]; j < G.ng[i]; ++j)
		{
                    B->gen[nr] = init_mat(dim,dim, "");
		    for (k = 0; k < dim; ++k)
		      for (l = 0; l < dim; ++l)
                          B->gen[nr]->array.SZ[k][l] = G.g[i][j][l][k];
		     ++nr;
		}
	}
        return(B);
}

/*********************************************************************\
|   writes the prime power decomposition
|   of G.ord[flags.STAB] *...* G.ord[G.dim-1] to B->divisors
\*********************************************************************/
void 
putord (group G, flagstruct flags, bravais_TYP *B)
{
	int	i, j, dim,  fac;

	dim = G.dim;
	for (i = 0; i < 100; ++i)
           B->divisors[i] = 0;
	B->divisors[1] = 1;
        B->order = 1;
	for (i = flags.STAB; i < dim; ++i)
	{
		fac = G.ord[i];
                B->order *= fac;
		for (j = 2; j <= dim+1  &&  fac != 1; ++j)
		{
			if (fac % j == 0)
			{
				++(B->divisors[j]);
				fac = fac / j;
				--j;
			}
		}
	}	
}

/*******************************************************\
| prints an isometry onto a matrix
\*******************************************************/
matrix_TYP *
putiso (int **X, int dim)
{
	int	i, j;
	matrix_TYP *M;

        M = init_mat(dim, dim, "");
	for (i = 0; i < dim; ++i)
		for (j = 0; j < dim; ++j)
                   M->array.SZ[i][j] = X[i][j];
        return(M);
}
