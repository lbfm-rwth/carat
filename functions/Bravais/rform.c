#include"typedef.h"
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

static double	***g, **x, **xt, **t, **t1, diff, diff1;
static int	dim, num;
static double	**sxt, **st, **st1, err, sdiff, sdiff1;




static int ggt(a, b)
int	a, b;
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



static void rmatfac(A, m, B, n)
double	**A, **B, m;
int	n;
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



static int rmatone(A, n)
double	**A;
int	n;
{
	int	i, j;


	for (i = 0; i < n; ++i)
	{
		if (fabs(A[i][i] - 1.0) > EPS)
			return 0;
		for (j = 0; j < n; ++j)
		{
			if (j != i  &&  fabs(A[i][j]) > EPS)
				return 0;
		}
	}
	return 1;
}


static void rmatadd(A, B, C, n)
double	**A, **B, **C;
int	n;
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



static void rmatmul(A, B, C, n)
double	**A, **B, **C;
int	n;
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


static void rmattrans(A, B, C, n)
double	**A, **B, **C;
int	n;
{
	double	**D, *d, c;
	int	i, j, k;


	if ((D = (double**)malloc(n * sizeof(double))) == 0)
		exit (2);
	for (i = 0; i < n; ++i)
	{
		if ((D[i] = (double*)malloc(n * sizeof(double))) == 0)
			exit (2);
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





static int gcd(f, dim)
int	**f, dim;
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

static int rond(X)
double	X;
{
	if (X >= 0)
		return (int)(2*fabs(X) + 1.0) / 2;
	else
		return -((int)(2*fabs(X) + 1.0) / 2);
}

static int chain(X, e)
double	X, e;
{
	double	*r, rest;
	int	i, *a, step, temp, nm = 0, dn = 1;
	double cond;

        /* inserted tilman 05/06/97 */
        if (e<EPS) e = EPS;
	cond = e;

	if ((r = (double*)malloc(sizeof(double))) == 0)
		exit (2);
	if ((a = (int*)malloc(sizeof(int))) == 0)
		exit (2);
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
		if ((r = (double*)realloc(r, (step+1) * sizeof(double))) == 0)
			exit (2);
		if ((a = (int*)realloc(a, step * sizeof(int))) == 0)
			exit (2);

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



static int dtoi(y, f, ep, dim)
double	**y, ep;
int	**f, dim;
{
	int	i, j, den, max;

	max = 1;
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			den = 1;
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
					den = chain(y[i][j], ep);
					if (den > max)
						max = den;
				}
			}
		}
	}
	return max;
}



static int iterate(step, x, g, num, dim)
double	**x, ***g;
int	step, num, dim;
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

	sdiff1 = 0;
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

	sdiff = 0;
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


static int symel(x, g, num, dim, eps, f)
double	**x, ***g;
int	num, dim, eps, **f;
{
	int	gc, i, j, k, den, maxden, step, neweps;

	if ((sxt = (double**)malloc(dim * sizeof(double*))) == 0)
		exit (2);
	if ((st = (double**)malloc(dim * sizeof(double*))) == 0)
		exit (2);
	if ((st1 = (double**)malloc(dim * sizeof(double*))) == 0)
			exit (2);
	for (i = 0; i < dim; ++i)
	{
		if ((sxt[i] = (double*)malloc(dim * sizeof(double))) == 0)
			exit (2);
		if ((st[i] = (double*)malloc(dim * sizeof(double))) == 0)
			exit (2);
		if ((st1[i] = (double*)malloc(dim * sizeof(double))) == 0)
			exit (2);
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



static double try()
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

	diff1 = 0;
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

	diff = 0;
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

matrix_TYP *rform(B, Banz, Fo, epsilon)
matrix_TYP **B;
int Banz;
matrix_TYP *Fo;
int epsilon;
{
	double	try();
	double	**save, **pres, *fac, factor;
	int	eps, best, d, i, j, k;
        matrix_TYP *erg;
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
             printf("error in rform: different dimesion of group elements\n");
             exit(3);
           }
        }
	num = Banz;
        erg = init_mat(dim,dim, "k");
	if ((g = (double***)malloc(num * sizeof(double**))) == 0)
		exit (2);
	if ((fac = (double*)malloc(num * sizeof(double))) == 0)
		exit (2);
	for (i = 0; i < num; ++i)
	{
		if ((g[i] = (double**)malloc(dim * sizeof(double*))) == 0)
			exit (2);
		for (j = 0; j < dim; ++j)
		{
			if ((g[i][j] = (double*)malloc(dim*sizeof(double)))==0)
				exit (2);
			for (k = 0; k < dim; ++k)
                          g[i][j][k] = (double)B[i]->array.SZ[k][j];
		}
	}

	if ((x = (double**)malloc(dim * sizeof(double*))) == 0)
		exit (2);
	if ((xt = (double**)malloc(dim * sizeof(double*))) == 0)
			exit (2);
	if ((t = (double**)malloc(dim * sizeof(double*))) == 0)
			exit (2);
	if ((t1 = (double**)malloc(dim * sizeof(double*))) == 0)
			exit (2);
	if ((save = (double**)malloc(dim * sizeof(double*))) == 0)
		exit (2);
	if ((pres = (double**)malloc(dim * sizeof(double*))) == 0)
		exit (2);
	for (i = 0; i < dim; ++i)
	{
		if ((x[i] = (double*)malloc(dim * sizeof(double))) == 0)
			exit (2);
		if ((xt[i] = (double*)malloc(dim * sizeof(double))) == 0)
			exit (2);
		if ((t[i] = (double*)malloc(dim * sizeof(double))) == 0)
			exit (2);
		if ((t1[i] = (double*)malloc(dim * sizeof(double))) == 0)
			exit (2);
		if ((save[i] = (double*)malloc(dim * sizeof(double))) == 0)
			exit (2);
		if ((pres[i] = (double*)malloc(dim * sizeof(double))) == 0)
			exit (2);
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
  	factor = try();
	for (i = 0; i < num; ++i)
	{
		for (j = 0; j < num; ++j)
		{
			if (j != i)
			{
				rmatfac(g[i], 1.0, pres, dim);
				rmatmul(g[i], g[j], t1, dim);
				rmatfac(t1, 1.0, g[i], dim);
				fac[j] = try();
				rmatfac(pres, 1.0, g[i], dim);
				rmatfac(save, 1.0, x, dim);
			}
		}
		fac[i] = factor;
		best = i;
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
