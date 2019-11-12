#include "typedef.h"
#include "utils.h"
#include"matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  rform.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

#define	DEFEPS		100
#define	DEFPRIME	5
#define	MAXDEN	1000000
#define	MAXSTEP	10000
#define	JUMP	5
#define	EPS	0.000001


static int 
ggt (int a, int b)
{
	int	r;

	if (b == 0)
		return(abs(a));

	while ((r = a % b) != 0)
	{
		a = b;
		b = r;
	}
	return (abs(b));
}



static void 
rmatfac (double **A, double m, double **B, int n)
{
	double	*a, *b;
	int	i, j;


	for (i = 0; i < n; ++i)
	{
		a = A[i];
		b = B[i];
		for (j = 0; j < n; ++j)
			b[j] = a[j] * m;
	}
}


static void 
rmatadd (double **A, double **B, double **C, int n)
{
	double	*a, *b, *c;
	int	i, j;


	for (i = 0; i < n; ++i)
	{
		a = A[i];
		b = B[i];
		c = C[i];
		for (j = 0; j < n; ++j)
			c[j] = a[j] + b[j];
	}
}



static void 
rmatmul (double **A, double **B, double **C, int n)
{
	double	*a, c;
	int	i, j, k;


	for (i = 0; i < n; ++i)
	{
		a = A[i];
		for (j = 0; j < n; ++j)
		{
			c = 0;
			for (k = 0; k < n; ++k)
			{
				if (a[k] != 0.0  &&  B[k][j] != 0.0)
					c += a[k] * B[k][j];
			}
			C[i][j] = c;
		}
	}
}


static void 
rmattrans (double **A, double **B, double **C, int n)
{
	double	**D, *d, c;
	int	i, j, k;


	D = (double**)xmalloc(n * sizeof(double *));
	for (i = 0; i < n; ++i)
	{
		D[i] = (double*)xmalloc(n * sizeof(double));
	}

	rmatmul(B, A, D, n);

	for (i = 0; i < n; ++i)
	{
		d = D[i];
		for (j = 0; j < n; ++j)
		{
			c = 0;
			for (k = 0; k < n; ++k)
				c += d[k] * B[j][k];
			C[i][j] = c;
		}
	}

	for (i = 0; i < n; ++i)
		free(D[i]);
	free(D);
}





static int 
gcd (int **f, int dim)
{
	int	n, i, j;

	n = f[0][0];
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			if (f[i][j] != 0)
				n = ggt(n, f[i][j]);
			if (n == 1)
				return 1;
		}
	}
	if (n == 0)
		return 1;
	else
		return n;
}

static int 
rond (double X)
{
	if (X >= 0)
		return (int)(2*fabs(X) + 1.0) / 2;
	else
		return -((int)(2*fabs(X) + 1.0) / 2);
}

static int 
chain (double X, double e)
{
	double	*r, rest;
	int	i, *a, step, temp, nm = 0, dn = 1;
	double cond;

        /* inserted tilman 05/06/97 */
        if (e<EPS) e = EPS;
	cond = e;

	r = (double*)xmalloc(sizeof(double));
	a = (int*)xmalloc(sizeof(int));
	temp = fabs(X);
	rest = fabs(X) - temp;
	step = 0;
	r[0] = rest;

	/* added the cond restriction because of numerical problems.
	for the formula of it see the expansion of 1/(x + epsilon) */
	while ((fabs(rest - (double)nm / (double)dn) > e ) &&
	       ((cond = fabs(cond / r[step]/ r[step]))< 1.0))
	{
		++step;
		r = (double*)xrealloc(r, (step+1) * sizeof(double));
		a = (int*)xrealloc(a, step * sizeof(int));

		a[step-1] = rond(1/r[step-1]);
		r[step] = 1/r[step-1] - a[step-1];
		nm = 1;
		dn = a[step-1];
		for (i = step-2; i >= 0; --i)
		{
			temp = dn;
			dn = a[i] * dn + nm;
			nm = temp;
		}

	}
	free(r);
	free(a);
	return abs(dn);
}



static int 
dtoi (double **y, int **f, double ep, int dim)
{
	int	i, j, max;

	max = 1;
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			/* changed 05/06/97 tilman from
			f[i][j] = y[i][j]; to: */
			f[i][j] = floor(y[i][j] + 0.5);
			if (fabs(y[i][j] - (double)f[i][j]) >= ep)
			{
				if (y[i][j] > 0  && 
						fabs(y[i][j]-(double)f[i][j]-1.0) < ep)
					f[i][j] += 1;
				else if (y[i][j] < 0  &&
						fabs(y[i][j]-(double)f[i][j]+1.0) < ep)
					f[i][j] -= 1;
				else
				{
					int den = chain(y[i][j], ep);
					if (den > max)
						max = den;
				}
			}
		}
	}
	return max;
}


static double	**xt, **t, **t1;
static double	**sxt, **st, **st1, err;

static int 
iterate (int step, double **x, double ***g, int num, int dim)
{
	double	factor;
	int	i, j, k;

	rmatfac(x, 1.0, st, dim);
	rmatfac(x, 1.0, sxt, dim);
	for (k = 0; k < JUMP; ++k)
	{
		++step;
		for (i = 0; i < num; ++i)
		{
			rmattrans(x, g[i], st1, dim);
			rmatadd(st, st1, st, dim);
		}
		for (i = 0; i < dim; ++i)
		{
			for (j = 0; j < dim; ++j)
			{
				st[i][j] *= 1/(double)(num+1);
				x[i][j] = st[i][j];
			}
		}
	}

	double sdiff1 = 0;
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			if (sdiff1 < fabs(sxt[i][j] - x[i][j]))
				sdiff1 = fabs(sxt[i][j] - x[i][j]);
		}
	}

	rmatfac(x, 1.0, sxt, dim);
	for (k = 0; k < JUMP; ++k)
	{
		++step;
		for (i = 0; i < num; ++i)
		{
			rmattrans(x, g[i], st1, dim);
			rmatadd(st, st1, st, dim);
		}
		for (i = 0; i < dim; ++i)
		{
			for (j = 0; j < dim; ++j)
			{
				st[i][j] *= 1/(double)(num+1);
				x[i][j] = st[i][j];
			}
		}
	}

	double sdiff = 0;
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			if (sdiff < fabs(sxt[i][j] - x[i][j]))
				sdiff = fabs(sxt[i][j] - x[i][j]);
		}
	}

	if (sdiff == 0 || sdiff1 == 0)
		err = 0.0;
	else
	{
		factor = pow(sdiff/sdiff1, 1.0/(double)JUMP);
		if (factor >= 1)
			err = 1;
		else
			err = sdiff*sdiff/(sdiff1-sdiff);
	/*	printf("step %d: sdiff = %.12lf  err: %.12lf  factor: %.12lf\n",
			step, sdiff, err, 1/factor);	*/
	}
	return step;
}


static int 
symel (double **x, double ***g, int num, int dim, int eps, int **f)
{
	int	gc, i, j, k, den, maxden, step, neweps;

	sxt = (double**)xmalloc(dim * sizeof(double*));
	st = (double**)xmalloc(dim * sizeof(double*));
	st1 = (double**)xmalloc(dim * sizeof(double*));
	for (i = 0; i < dim; ++i)
	{
		sxt[i] = (double*)xmalloc(dim * sizeof(double));
		st[i] = (double*)xmalloc(dim * sizeof(double));
		st1[i] = (double*)xmalloc(dim * sizeof(double));
	}
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
			st[i][j] = x[i][j];
	}

	step = 0;
	err = 1.0;
	step = iterate(step, x, g, num, dim);
	while (err > 0.3/(double)eps  &&  step < MAXSTEP)
		step = iterate(step, x, g, num, dim);
	if (step >= MAXSTEP)
	{
		printf("Verfahren hat nach %d Schritten nicht konvergiert !\n", step);
		exit(3);
		return 0;
	}
	den = 0;
	while (den != 1)
	{
		den = dtoi(x, f, 2*err, dim);
		if (den != 1)
		{
			maxden = den+1;
			while (maxden > den  &&  maxden < MAXDEN)
			{
				neweps = (maxden < eps) ? eps*maxden : maxden*maxden;

				/* inserted tilman 09/06/97 */
				if (1/(double) neweps < EPS) neweps = 1/EPS;

				while (err > 0.3/(double)neweps  &&  step < MAXSTEP)
					step = iterate(step, x, g, num, dim);
				if (step >= MAXSTEP)
				{
					printf("Verfahren hat nach %d Schritten nicht konvergiert !\n", step);
					exit(3);
					return 0;
				}
				den = maxden;
			/*	printf("Versuch mit %d\n", den);	*/
				maxden = dtoi(x, f, 2*err, dim);
			}
			rmatfac(x, (double)maxden, x, dim);
		/*	printf("Nenner mit %d multipliziert\n", maxden);*/
			err *= maxden;
		}
	}

	gc = gcd(f, dim);
/*	printf("ggt %d\n", gc);
*/	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			f[i][j] /= gc;
			st[i][j] = (double)f[i][j];
		}
	}
	for (i = 0; i < num; ++i)
	{
		rmattrans(st, g[i], st1, dim);
		for (j = 0; j < dim; ++j)
		{
			for (k = 0; k < dim; ++k)
			{
				if (st1[j][k] != st[j][k])
				{
					printf("Matrix tuts nicht !\n");
					exit(3);
					return 0;
				}
			}
		}
	}
  for(i=0;i<dim;i++)
    { free(sxt[i]); free(st[i]); free(st1[i]);}
  free(sxt); free(st); free(st1);
	return (step);
}


/*
static double 
try (double **x, double ***g, const int num, const int dim)
{
	double	conv;
	int	i, j, k;

	rmatfac(x, 1.0, t, dim);
	for (k = 0; k < 4; ++k)
	{
		for (i = 0; i < num; ++i)
		{
			rmattrans(x, g[i], t1, dim);
			rmatadd(t, t1, t, dim);
		}
		for (i = 0; i < dim; ++i)
		{
			for (j = 0; j < dim; ++j)
			{
				t[i][j] *= 1/(double)(num+1);
				x[i][j] = t[i][j];
			}
		}
	}

	rmatfac(x, 1.0, xt, dim);
	for (k = 0; k < 3; ++k)
	{
		for (i = 0; i < num; ++i)
		{
			rmattrans(x, g[i], t1, dim);
			rmatadd(t, t1, t, dim);
		}
		for (i = 0; i < dim; ++i)
		{
			for (j = 0; j < dim; ++j)
			{
				t[i][j] *= 1/(double)(num+1);
				x[i][j] = t[i][j];
			}
		}
	}

	double diff1 = 0;
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			if (diff1 < fabs(xt[i][j] - x[i][j]))
				diff1 = fabs(xt[i][j] - x[i][j]);
		}
	}

	rmatfac(x, 1.0, xt, dim);
	for (k = 0; k < 3; ++k)
	{
		for (i = 0; i < num; ++i)
		{
			rmattrans(x, g[i], t1, dim);
			rmatadd(t, t1, t, dim);
		}
		for (i = 0; i < dim; ++i)
		{
			for (j = 0; j < dim; ++j)
			{
				t[i][j] *= 1/(double)(num+1);
				x[i][j] = t[i][j];
			}
		}
	}

	double diff = 0;
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			if (diff < fabs(xt[i][j] - x[i][j]))
				diff = fabs(xt[i][j] - x[i][j]);
		}
	}

	conv = pow(diff/diff1, 1.0/3.0);
	return (1/conv);
}
*/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *rform(B, Banz, Fo, epsilon)
@ matrix_TYP **B;
@ int Banz;
@ matrix_TYP *Fo;
@ int epsilon;
@
@ If G denotes the group generated by B[0],...,B[Banz-1],
@ then 'rform' calculates a matrix  A that is the sum over all g in G of
@   g^{tr} * Fo * g
@ This matrix is divided by the gcd of its entries and returned.
@ The algorithm uses an appoximation and the precession is 1/epsilon.
@ If epsilon < 100, the function sets epsilon = 100.
@---------------------------------------------------------------------------
@
\**************************************************************************/

matrix_TYP *
rform (matrix_TYP **B, int Banz, matrix_TYP *Fo, int epsilon)
{
	double	**save, **pres, *fac;
	double  ***g, **x;
	int	eps, d, i, j, k;
        matrix_TYP *erg;
    int num, dim;
    int test;

	eps = epsilon;
	if (eps < DEFEPS)
		eps = DEFEPS;
	dim = B[0]->cols;
        if(B[0]->rows != dim)
        {
          printf("error in rform: non-square matrix in group 'B'\n");
          exit(3);
        }
        for(i=1;i<Banz;i++)
        {
           if(B[i]->rows != dim || B[i]->cols != dim)
           {
             printf("error in rform: different dimension of group elements\n");
             exit(3);
           }
        }
	num = Banz;
        erg = init_mat(dim,dim, "k");
	g = (double***)xmalloc(num * sizeof(double**));
	fac = (double*)xmalloc(num * sizeof(double));
	for (i = 0; i < num; ++i)
	{
		g[i] = (double**)xmalloc(dim * sizeof(double*));
		for (j = 0; j < dim; ++j)
		{
			g[i][j] = (double*)xmalloc(dim*sizeof(double));
			for (k = 0; k < dim; ++k)
                          g[i][j][k] = (double)B[i]->array.SZ[k][j];
		}
	}

	x = (double**)xmalloc(dim * sizeof(double*));
	xt = (double**)xmalloc(dim * sizeof(double*));
	t = (double**)xmalloc(dim * sizeof(double*));
	t1 = (double**)xmalloc(dim * sizeof(double*));
	save = (double**)xmalloc(dim * sizeof(double*));
	pres = (double**)xmalloc(dim * sizeof(double*));
	for (i = 0; i < dim; ++i)
	{
		x[i] = (double*)xmalloc(dim * sizeof(double));
		xt[i] = (double*)xmalloc(dim * sizeof(double));
		t[i] = (double*)xmalloc(dim * sizeof(double));
		t1[i] = (double*)xmalloc(dim * sizeof(double));
		save[i] = (double*)xmalloc(dim * sizeof(double));
		pres[i] = (double*)xmalloc(dim * sizeof(double));
	}
        d = Fo->cols;
	if (d != dim)
	{
		printf("Fehler : Unkompatible Dimensionen\n");
		exit (3);
	}
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
			x[i][j] = (double) Fo->array.SZ[i][j];
	}
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
			save[i][j] = x[i][j];
	}

        /*
  	int factor = try(x, g, num, dim);
	for (i = 0; i < num; ++i)
	{
		for (j = 0; j < num; ++j)
		{
			if (j != i)
			{
				rmatfac(g[i], 1.0, pres, dim);
				rmatmul(g[i], g[j], t1, dim);
				rmatfac(t1, 1.0, g[i], dim);
				fac[j] = try(x, g, num, dim);
				rmatfac(pres, 1.0, g[i], dim);
				rmatfac(save, 1.0, x, dim);
			}
		}
		fac[i] = factor;
		int best = i;
		for (j = 0; j < num; ++j)
		{
			if (fac[j] > fac[best])
				best = j;
		}
		if (best != i)
		{
			rmatmul(g[i], g[best], t1, dim);
			rmatfac(t1, 1.0, g[i], dim);
			i = -1;
		}
		factor = fac[best];
	}
	*/ 

	 test = symel(x, g, num, dim, eps, erg->array.SZ);

   for(i=0;i<num;i++)
   {
     for(j=0;j<dim;j++)
       free(g[i][j]);
     free(g[i]);
   }
   free(g);
   for(i=0;i<dim;i++)
   { 
     free(x[i]); free(xt[i]); free(t[i]); free(t1[i]);
     free(save[i]); free(pres[i]);
   }
   free(x); free(xt); free(t); free(t1);
   free(save); free(pres);
   free(fac);

  if(test != 0)
    return(erg);
  else
  { free_mat(erg); return(NULL);}
}
