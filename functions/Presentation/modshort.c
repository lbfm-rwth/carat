
/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: modshort.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

/**********************************************
 Programm berechnet die Vektoren in vi+T bis zur Laenge d
 rationaler Startvektor moeglich 
   in mat steht die Grammatrix TAT_tr des Gitters, (bzgl eines invarianten 
   Skalarproduktes A der Punktgruppe -- option: eingeben oder mit RFORM --),
   restvec ist der Vektor um den das Gitter verschoben sein soll. 
   ist restvector rational mit Nenner kgv so ist grn=kgv^2*mat , suche kuerzeste
   Vektoren in kgv*restvec+kgv*T und dividiere das Ergebnis durch kgv
	matrix_TYP *neumodshort_vectors(mat, restvector, d)
**********************************************/
#include"typedef.h"
#include"matrix.h"

#define	EPS	0.001

/*static	int	n, max, *vec, anzahl,kgv;*/
static	int	n, max, *vec, anzahl,kgv;
static	int	**grn;		/* Gram matrix */
static	int	**ba;		/* basis transformation matrix */
static	double	**mo, *ge;
static matrix_TYP *SV;
static int SV_size, SV_ext;
static int *restvec;
static double *restvek, rk;
static int *rv;
static double scapr(int i, int j);
static double orth(int i);
static void modellmachen(int a, int e);
static int round( double r);
static void red(int k, int l);
static int interchange(int k);
static void reduzieren(int k);
static void vecschr(int m, double d);
static void shrt(int c, double damage);

void anhalter()
{
  int i;
  i = 1;
}

static double scapr(i, j)
int	i, j;
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
                anhalter();
		exit(3);
	}
        else
	{
	r = r / ge[j];
	}
	return r;
}

static double orth(i)
int	i;
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


static void modellmachen(a, e)
int	a, e;
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

static int round(r)
double	r;
{
	int	i;

	if (r >= 0)
		i = (int)(2*r + 1) / 2;
	else
		i = -(int)(2*(-r) + 1) / 2;
	return i;
}
	
static void red(k, l)
int	k, l;
{
	double	r, *mok, *mol;
	int	i, ir, *bak, *bal, *grnk, *grnl;

	r = mo[k][l];
	if (2*r > 1.0  ||  2*r < -1.0)
	{
		ir = round(r);
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

static int interchange(k)
int	k;
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

static void reduzieren(k)
int	k;
{
	int	l;

	for (l = k-2; l >= 0; --l)
		red(k, l);
}


/* prints the vector corresponding to vec */
/* in the model and its length            */

static void vecschr(m, d)
int m;
double	d;
{
	int	i, j;
        int l;
        int *v;

        l = round(d);
        if(l<max)             /*(hier wird nur der kuerzeste Vektor bestimmt)*/
        {
         max = l;
         for(i=0;i<SV->rows;i++)
           free(SV->array.SZ[i]);
         SV_size = SV_ext;
         SV->rows = 0;
         if((SV->array.SZ = (int **)realloc(SV->array.SZ, SV_size *sizeof(int *))) == 0)
                  exit(2);
	 anzahl = 0;
        }
	 anzahl++;
        if(SV->rows == SV_size)
        {
           SV_size += SV_ext;
           if((SV->array.SZ = (int **)realloc(SV->array.SZ, SV_size *sizeof(int *))) == 0)
                exit(2);
        }
        if((v = (int *)malloc((n+1) *sizeof(int))) == 0)
            exit(2);
        for (i = 0; i < m; ++i)
        {
		  v[i] = 0;
		  for (j = 0; j < m; ++j)
			  v[i] += vec[j] * ba[j][i];
        }
        v[m] = l;
/*
for (i = 0; i < m; ++i)
printf("%d\n",v[i]);
printf("hallo\n");
*/
        SV->array.SZ[SV->rows] = v;
        SV->rows++;
 /*       SV->kgv = kgv;	dividiere das Ergebnis durch restvector->kgv*/
}



 /* recursion for finding shortest vectors */

static void shrt(c, damage)
int	c;
double	damage;
{
	double	x, gec;
	int	i, j;
/*
printf("c= %d\n",c);
printf("damage = %f\n", damage);
*/
	if (c == -1)
	{
			vecschr(n, damage);
	}
	else
	{
		x = 0.0;
		for (j = c+1; j < n; ++j)
			x += ((double)vec[j]-restvek[j]) * mo[j][c];
       		x -= restvek[c];
		i = -round(x);	
		gec = ge[c];
/*
printf("gec=%f\n",gec);
*/
  	        if (gec*(x+i)*(x+i) + damage < max + EPS)
		{
			while ((gec*(x+i)*(x+i) + damage <= max + EPS) || 
				(x + i <= 0))
				++i;
			--i;
			while ((gec*(x+i)*(x+i) + damage < max + EPS))
			{
 		 		vec[c] = i;
		 		shrt(c-1, gec*(x+i)*(x+i)+damage);
				--i;
			}
		}
	}
}


/* in mat steht die Grammatrix TAT_tr des Gitters, (bzgl eines invarianten 
   Skalarproduktes A der Punktgruppe -- option: eingeben oder mit RFORM --),
   restvec ist der Vektor um den das Gitter verschoben sein soll. 
   ist restvector rational mit Nenner k so ist grn=k^2*mat , suche kuerzeste
   Vektoren in k*restvec+k*T und dividiere das Ergebnis durch k*/

/**************************************************************************\
@--------------------------------------------------------------------------
@ matrix_TYP *modshort_vectors(mat, restvector,d)
@
@ matrix_TYP *mat:         grammatrix of a lattice L
@ matrix_TYP *restvector:  a (possibly rational) vector
@ int d:                   given length
@
@ calculates the vectors in restvector+L up to length d and returns these as
@ the rows of a matrix.
@
@--------------------------------------------------------------------------
\**************************************************************************/
matrix_TYP *modshort_vectors(mat, restvector,d)
matrix_TYP *mat;
matrix_TYP *restvector;
int d;
{
	int	i, j, ak, bk;
	double	cst = 0.75;
        matrix_TYP *erg;

	anzahl = 0;
        n = mat->rows;
        kgv = restvector->kgv;
        restvec = (int *)malloc(n *sizeof(int));
        for(i = 0; i < n; i++)
         restvec[i] = restvector->array.SZ[0][i]; 

        SV_ext = n * (n+1)/2;

        SV = init_mat(1, n+1, "");
        free(SV->array.SZ[0]);
        SV->cols--;
        SV->rows = 0;
        
          if((SV->array.SZ = (int **)realloc(SV->array.SZ, SV_ext *sizeof(int *))) == 0)
                exit(2);
          SV_size = SV_ext;
        
	if ((grn = (int**)malloc(n * sizeof(int*))) == 0)
		exit (2);
	if ((ba = (int**)malloc(n * sizeof(int*))) == 0)
		exit (2);
	if ((mo = (double**)malloc(n * sizeof(double*))) == 0)
		exit (2);
	if ((ge = (double*)malloc(n * sizeof(double))) == 0)
		exit (2);
	if ((vec = (int*)malloc(n * sizeof(int))) == 0)
		exit (2);
	if ((rv = (int*)malloc(n * sizeof(int))) == 0)
		exit (2);
	if ((restvek = (double*)malloc(n * sizeof(double))) == 0)
		exit (2);
	for (i = 0; i < n; ++i)
	{
		vec[i] = 0;
		rv[i] = restvec[i];
		if ((grn[i] = (int*)malloc(n * sizeof(int))) == 0)
			exit (2);
		if ((ba[i] = (int*)malloc(n * sizeof(int))) == 0)
			exit (2);
		if ((mo[i] = (double*)malloc(n * sizeof(double))) == 0)
			exit (2);
		for (j = 0; j < n; ++j)
		{
			ba[i][j] = 0;
			mo[i][j] = 0.0;
		}
		ba[i][i] = 1;
		mo[i][i] = 1.0;

/* read the Gram matrix (can be square or lower triangular) */
          	for(j = 0; j <= i; j++)
		{
            		grn[i][j] = (mat->array.SZ[i][j]);
			grn[j][i] = grn[i][j];
		}
	}

	rk = ((double)kgv);
	max = d; 
/* vgl... min = 0;  */
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
/*
for(i=0; i<n; i++)
printf("ge = %f ",ge[i]);
printf("\n ");
*/        
	for(i=0; i<n; i++)
	restvek[i] = ((double)rv[i])/rk;

	shrt(n-1, 0.0);	/* recursively calculate the vectors up to length max */

     for(i=0;i<n;i++)
     {
         free(grn[i]);
         free(ba[i]);
         free(mo[i]);
     }
     free(grn); free(ba); free(mo); free(ge); free(vec);
     free(rv); free(restvek);
     free(restvec);/*4.3.anne*/
     
     Check_mat(SV);
     erg = SV;
     SV = NULL;
     return(erg);
}/**modshort_vectors()**/
