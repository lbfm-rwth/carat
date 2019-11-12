#include "typedef.h"
#include "utils.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: shortest.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

#define	EPS	0.001

static	int	n, max, *vec, anzahl, con;
static	int	**grn;		/* Gram matrix */
static	int	**ba;		/* basis transformation matrix */
static	double	**mo, *ge;
static matrix_TYP *SV;
static int SV_size, SV_ext;


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
        if(l<max)
        {
         max = l;
         for(i=0;i<SV->rows;i++)
           free(SV->array.SZ[i]);
         SV_size = SV_ext;
         SV->rows = 0;
         SV->array.SZ = (int **)xrealloc(SV->array.SZ, SV_size *sizeof(int *));
          anzahl = 0;
        }
        anzahl++;
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


 /* recursion for finding shortest vectors */

static void 
shrt (int c, double damage)
{
	double	x, gec;
	int	i, j;

	if (c == -1)
	{
		for (i = 0; i < n  &&  vec[i] == 0; ++i);
		if (i == n)
			con = 1;
		else
		{
			vecschr(n, damage);
		}
	}
	else
	{
		x = 0.0;
		for (j = c+1; j < n; ++j)
			x += vec[j] * mo[j][c];
		i = -iround(x);	
		gec = ge[c];
		if (gec*(x+i)*(x+i) + damage < max + EPS)
		{
			while ((gec*(x+i)*(x+i) + damage <= max + EPS) ||
				(x + i <= 0))
				++i;
			--i;
			while ((gec*(x+i)*(x+i) + damage < max + EPS) && 
				con == 0)
			{
				vec[c] = i;
				shrt(c-1, gec*(x+i)*(x+i) + damage);
				--i;
			}
		}
	}
}




/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *shortest(mat, min_norm)
@ matrix_TYP *mat;
@ int  *min_norm;
@ 
@ The function callulates the shortest vectors of the integral, symmetric,
@ positiv definite matrix 'mat'.
@ The minmum
@ 
@    m(mat) := min{  v * mat *v^{tr} | v in Z^n and v not 0}
@ 
@ is returned in the pointer 'min_norm'. The shortest vectors are the vectors
@ with n(x) = min(x).
@ They are stored as rows in the returned matrix called 'SV' in the following.
@ For 'SV' n+1 columns are allocated. The last column contains the norm n(x).
@ Although n+1 columns are allocated, SV->cols = n.
@ 
@ Warning: if 'mat' is not positiv definite, the function runs into an 
@         infinite loop.
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
shortest (matrix_TYP *mat, int *min_norm)
{
	int	i, j, ak, bk;
	double	cst = 0.75;
        matrix_TYP *erg;

        anzahl = 0;
        con = 0;
        n = mat->rows;
        SV_ext = n * (n+1)/2;
        SV = init_mat(1, n+1, "");
        SV->cols--;
        SV->rows = 0;
        free(SV->array.SZ[0]);
        /* inserted the next line 15/1/97 tilman */
        free(SV->array.SZ);

        SV->array.SZ = (int **)xmalloc(SV_ext *sizeof(int *));
        SV_size = SV_ext;
	grn = (int**)xmalloc(n * sizeof(int*));
	ba = (int**)xmalloc(n * sizeof(int*));
	mo = (double**)xmalloc(n * sizeof(double*));
	ge = (double*)xmalloc(n * sizeof(double));
	vec = (int*)xmalloc(n * sizeof(int));
	for (i = 0; i < n; ++i)
	{
		vec[i] = 0;
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
        max = mat->array.SZ[0][0];
        for(i=1;i<n;i++)
        {
            if(mat->array.SZ[i][i] < max)
              max = mat->array.SZ[i][i];
        }

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

	shrt(n-1, 0.0);	/* recursively calculate the vectors up to length max */

     for(i=0;i<n;i++)
     {
         free(grn[i]);
         free(ba[i]);
         free(mo[i]);
     }
     free(grn); free(ba); free(mo); free(ge); free(vec);

     erg = SV;
     SV = NULL;
     *min_norm = max;

     return(erg);
}
