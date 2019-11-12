#include "typedef.h"
#include "utils.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: rest_short.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

#define	EPS	0.001

static	int	n, min, *vec, anzahl;
static	int	**grn;		/* Gram matrix */
static	int	**ba;		/* basis transformation matrix */
static	double	**mo, *ge;
static matrix_TYP *SV;
static int SV_size, SV_ext;
static int break_opt, f_opt, c_opt;
static int *rv;
static double *rvd, rk, max;


static double 
scapr (int i, int j)
{
	double	r, *moi, *moj;
	int	l;

	r = grn[i][j];
	moi = mo[i];
	moj = mo[j];
	for (l = 0; l <= j-1; ++l)
		r -= ge[l] * moi[l] * moj[l];
	if (ge[j] == 0)
	{
		printf("Error: found norm 0 vector\n");
		exit(3);
	}
	r = r / ge[j];
	return r;
}

static double 
orth (int i)
{
	double	r, *moi;
	int	l;

	r = grn[i][i];
	moi = mo[i];
	for (l = 0; l <= i-1; ++l)
		r -= ge[l] * moi[l] * moi[l];
	return r;	
}

/*---------------------------------------------------------*\
|   makes model of the lattice (with same scalarproducts)
\*__________________________________________________________*/

static void 
modellmachen (int a, int e)
{
	int	i, j;

	if (a == 0)
	{
		ge[0] = grn[0][0];
		i = 1;
	}
	else
		i = a;
	while (i <= e)
	{
		for (j = 0; j < i; ++j)
			mo[i][j] = scapr(i, j);
		ge[i] = orth(i);
		++i;
	}
}

static int 
iround (double r)
{
	int	i;

	if (r >= 0)
		i = (int)(2*r + 1) / 2;
	else
		i = -(int)(2*(-r) + 1) / 2;
	return i;
}
	
static void 
red (int k, int l)
{
	double	r, *mok, *mol;
	int	i, ir, *bak, *bal, *grnk, *grnl;

	r = mo[k][l];
	if (2*r > 1.0  ||  2*r < -1.0)
	{
		ir = iround(r);
		mok = mo[k];
		bak = ba[k];
		grnk = grn[k];
		mol = mo[l];
		bal = ba[l];
		grnl = grn[l];
		for (i = 0; i < n; ++i)
		{
			mok[i] -= ir * mol[i];
			bak[i] -= ir * bal[i];
			grnk[i] -= ir * grnl[i];	
		}
                rv[l] += ir * rv[k];
		for (i = 0; i < n; ++i)
			grn[i][k] -= ir * grn[i][l];	
	}
}

static int 
interchange (int k)
{
	int	i, z, *zp;

	zp = ba[k];
	ba[k] = ba[k-1];
	ba[k-1] = zp;

        z = rv[k];
        rv[k] = rv[k-1];
        rv[k-1] = z;

	zp = grn[k];
	grn[k] = grn[k-1];
	grn[k-1] = zp;

	for (i = 0; i < n; ++i)
	{
		z = grn[i][k];
		grn[i][k] = grn[i][k-1];
		grn[i][k-1] = z;
	}
	if (k > 1)
		return k-1;
	else
		return k;
}

static void 
reduzieren (int k)
{
	int	l;

	for (l = k-2; l >= 0; --l)
		red(k, l);
}


/* prints the vector corresponding to vec */
/* in the model and its length            */

static void 
vecschr (int m, double d)
{
	int	i, j;
        int l;
        int *v;

        l = iround(d);
        if(l >= min)
        {
          anzahl++;
          if(f_opt == TRUE)
            break_opt = TRUE;
          if(c_opt == TRUE)
             return;
          if(SV->rows == SV_size)
          {
             SV_size += SV_ext;
             SV->array.SZ = (int **)xrealloc(SV->array.SZ, SV_size *sizeof(int *));
          }
          v = (int *)xmalloc((n+1) *sizeof(int));
	  for (i = 0; i < m; ++i)
	  {
		  v[i] = 0;
		  for (j = 0; j < m; ++j)
			  v[i] += vec[j] * ba[j][i];
	  }
          v[m] = l;
          SV->array.SZ[SV->rows] = v;
          SV->rows++;
        }
}


 /* recursion for finding shortest vectors */

static void 
shrt (int c, double damage)
{
	double	x, gec;
	int	i, j;

	if (c == -1)
	{
			vecschr(n, damage);
                        if(break_opt == TRUE)
                          return;
	}
	else
	{
		x = 0.0;
		for (j = c+1; j < n; ++j)
			x += (vec[j] -rvd[j]) * mo[j][c];
                x -= rvd[c];
		i = -iround(x);	
		gec = ge[c];
		if (gec*(x+i)*(x+i) + damage < max + EPS)
		{
			while ((gec*(x+i)*(x+i) + damage <= max + EPS) ||
				(x + i <= 0))
				++i;
			--i;
			while ((gec*(x+i)*(x+i) + damage < max + EPS))
			{
				vec[c] = i;
                                if(break_opt == TRUE)
                                    return;
				shrt(c-1, gec*(x+i)*(x+i) + damage);
				--i;
			}
		}
	}
}






/**************************************************************************\
@---------------------------------------------------------------------------
matrix_TYP *short_vectors(mat, length, lengthmin, find_opt, count_opt, anz)
matrix_TYP *mat;
int length, lengthmin, find_opt, count_opt, *anz;

@ matrix_TYP *rest_short(mat, restvec, rkgv, zaehler, nenner,
@                        find_opt, count_opt, anz)
@ matrix_TYP *mat;
@ int *restvec, rkgv;
@ int zaehler, nenner, find_opt, count_opt, *anz;
@
@ The matrix mat must be a symmetric, positiv definite nxn-Matrix.
@ The function calculates all x in Z^n  with 
@ 
@          lengthmin <= n(x - 1/rkgv * restvec)  <= length
@ where
@          n(y) := y * mat * y^{tr}.
@ 
@ In the pointer 'anz' the number of all these vectors is returned.
@ The vectors x are the rows of the returned matrix, called 'SV' in the
@ following.
@ So 'SV' is an integral matrix with 'anz' rows and 'n' columns.
@ For 'SV' n+1 columns are allocated. The last column contains the norm n(x).
@ Although n+1 columns are allocated, SV->cols = n.
@ 
@ If 'findoption == 1' the Algorithm stops,
@ if a vector of desired norm is found.
@ In this case the Matrix 'SV' has only one row.
@ 
@ If 'count_option == 1',  shortest_vectors returns a zero-pointer,
@ only the number of found vectors is returned in the pointer 'anz'.
@ 
@ Warning: if 'mat' is not positiv definite, the function runs into an 
@          infinite loop.
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
rest_short (matrix_TYP *mat, int *restvec, int rkgv, int zaehler, int nenner, int find_opt, int count_opt, int *anz)
{
	int	i, j, ak, bk;
	double	cst = 0.75;
        matrix_TYP *erg;

        anzahl = 0;
        n = mat->rows;
        SV_ext = n * (n+1)/2;
        break_opt = FALSE;
        f_opt = find_opt;
        c_opt = count_opt;

        SV = init_mat(1, n+1, "");
        free(SV->array.SZ[0]);
        SV->cols--;
        SV->rows = 0;
        if(f_opt == TRUE)
          SV_ext = 1;
        else
        {
          SV->array.SZ = (int **)xrealloc(SV->array.SZ, SV_ext *sizeof(int *));
          SV_size = SV_ext;
        }
	grn = (int**)xmalloc(n * sizeof(int*));
	ba = (int**)xmalloc(n * sizeof(int*));
	mo = (double**)xmalloc(n * sizeof(double*));
	ge = (double*)xmalloc(n * sizeof(double));
	vec = (int*)xmalloc(n * sizeof(int));
	rv = (int*)xmalloc(n * sizeof(int));
	rvd = (double*)xmalloc(n * sizeof(double));
	for (i = 0; i < n; ++i)
	{
		vec[i] = 0;
                rv[i] = restvec[i];
		grn[i] = (int*)xmalloc(n * sizeof(int));
		ba[i] = (int*)xmalloc(n * sizeof(int));
		mo[i] = (double*)xmalloc(n * sizeof(double));
		for (j = 0; j < n; ++j)
		{
			ba[i][j] = 0;
			mo[i][j] = 0.0;
		}
		ba[i][i] = 1;
		mo[i][i] = 1.0;

/* read the Gram matrix (can be square or lower triangular) */
		for (j = 0; j <= i; ++j)
		{
			grn[i][j] = mat->array.SZ[i][j];
			grn[j][i] = grn[i][j];
		}
	}
        rk = ((double) rkgv);
	max = (double) zaehler;	/* read the maximal length of the vectors */
        max *= ((double) mat->kgv);
        if(nenner != 0)
        max = max /((double) nenner);
        min = 0;

	bk = 1;
	while (bk < n)
	{
		ak = bk - 1;
		modellmachen(ak, bk);
		red(bk, bk-1);
		if (ge[bk] < ((cst - mo[bk][bk-1]*mo[bk][bk-1]) * ge[bk-1]))
		{
			bk = interchange(bk);
			modellmachen(bk-1, bk);
			if (bk > 1)
				--bk;
		}
		else
		{
			if (bk > 1)
				reduzieren(bk);
			++bk;
		}
	}
        for(i=0;i<n;i++)
          rvd[i] = ((double) rv[i])/rk;

	shrt(n-1, 0.0);	/* recursively calculate the vectors up to length max */
        *anz = anzahl;

     for(i=0;i<n;i++)
     {
         free(grn[i]);
         free(ba[i]);
         free(mo[i]);
     }
     free(grn); free(ba); free(mo); free(ge); free(vec);
     free(rv); free(rvd);

     erg = SV;
     SV = NULL;
     return(erg);
}
