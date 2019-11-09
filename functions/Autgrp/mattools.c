/*****	This file includes some basic routines 
	for matrix and vector operations	*****/
#include "typedef.h"
#include "types.h"

/********************************************************************\
| multiplies 1xn-vector x with nxn-matrix A and puts result on y
\********************************************************************/
void 
vecmatmul (int *x, int **A, int n, int *y)
{
	int	i, j, xi, *Ai;

	for (j = 0; j < n; ++j)
		y[j] = 0;
	for (i = 0; i < n; ++i)
	{
		if ((xi = x[i]) != 0)
		{
			Ai = A[i];
			for (j = 0; j < n; ++j)
				y[j] += xi * Ai[j];
		}
	}
}

/********************************************************************\
|	multiplies nxn-matrices A and B and puts result	on C
\********************************************************************/
void 
matmul (int **A, int **B, int n, int **C)
{
	int	i, j, k, *Ai, *Bk, *Ci, Aik;

	for (i = 0; i < n; ++i)
	{
		Ai = A[i];
		Ci = C[i];
		for (j = 0; j < n; ++j)
			Ci[j] = 0;
		for (k = 0; k < n; ++k)
		{
			if ((Aik = Ai[k]) != 0)
			{
				Bk = B[k];
				for (j = 0; j < n; ++j)
					Ci[j] += Aik * Bk[j];
			}
		}
	}
}

/********************************************************************\
|	returns scalar product of 1xn-vectors x and y
|	with respect to Gram-matrix F
\********************************************************************/
int 
scp (int *x, int **F, int *y, int n)
{
	int	i, j, sum, xi, *Fi;

	sum = 0;
	for (i = 0; i < n; ++i)
	{
		if ((xi = x[i]) != 0)
		{
			Fi = F[i];
			for (j = 0; j < n; ++j)
				sum += xi * Fi[j] * y[j];
		}
	}
	return(sum);
}

/********************************************************************\
|  returns standard scalar product of 1xn-vectors x and y, heavily used
\********************************************************************/
int 
sscp (int *x, int *y, int n)
{
	int	i, sum;

	sum = 0;
	for (i = 0; i < n; ++i)
		sum += *(x++) * *(y++);
	/*	sum += x[i] * y[i];	*/
	return(sum);
}

/********************************************************************\
|	computes X = B*A^-1 modulo the prime p, for nxn-matrices A and B,
|	where A is invertible, A and B are modified !!!	
\********************************************************************/
void 
psolve (int **X, int **A, int **B, int n, int p)
{
	int	i, j, k, tmp, sum, ainv, *Xi, *Bi;

/* convert A to lower triangular matrix and change B accordingly */
	for (i = 0; i < n-1; ++i)
	{
		for (j = i; A[i][j] % p == 0; ++j);
		if (j == n)
		{
			fprintf(stderr, "Error: matrix is singular modulo %d\n", p);
			exit (3);
		}
		if (j != i)
/* interchange columns i and j such that A[i][i] is != 0 */
		{
			for (k = i; k < n; ++k)
			{
				tmp = A[k][i];
				A[k][i] = A[k][j];
				A[k][j] = tmp;
			}
			for (k = 0; k < n; ++k)
			{
				tmp = B[k][i];
				B[k][i] = B[k][j];
				B[k][j] = tmp;
			}
		}
		pgauss(i, A, B, n, p);
	}
/* compute X recursively */
	for (i = 0; i < n; ++i)
	{
		Xi = X[i];
		Bi = B[i];
		for (j = n-1; j >= 0; --j)
		{
			sum = 0;
			for (k = n-1; k > j; --k)
				sum = (sum + Xi[k] * A[k][j]) % p;
/* ainv is the inverse of A[j][j] modulo p */
			for (ainv = 1; abs(A[j][j]*ainv-1) % p != 0; ++ainv);
			Xi[j] = (Bi[j] - sum) * ainv % p;
/* make sure that -p/2 < X[i][j] <= p/2 */
			if (2*Xi[j] > p)
				Xi[j] -= p;
			else if (2*Xi[j] <= -p)
				Xi[j] += p;
		}
	}
}

/********************************************************************\
|	clears row nr. r of A assuming that the 
|	first r-1 rows of A are already cleared 
|	and that A[r][r] != 0 modulo p, 
|	A and B	are changed !!!	
\********************************************************************/
void 
pgauss (int r, int **A, int **B, int n, int p)
{
	int	i, j, f, ainv, *Ar;

	Ar = A[r];
/* ainv is the inverse of A[r][r] modulo p */
	for (ainv = 1; abs(Ar[r]*ainv-1) % p != 0; ++ainv);
	for (j = r+1; j < n; ++j)
	{
		if (Ar[j] % p != 0)
		{
			f = Ar[j] * ainv % p;
			for (i = r+1; i < n; ++i)
				A[i][j] = (A[i][j] - f * A[i][r]) % p;
			for (i = 0; i < n; ++i)
				B[i][j] = (B[i][j] - f * B[i][r]) % p;
			Ar[j] = 0;
		}
	}
}

/********************************************************************\
|	checks, whether n is a prime
\********************************************************************/
int 
isprime (int n)
{
	int	i;

	for (i = 2; i <= n/i; ++i)
	{
		if (n % i == 0)
			return 0;
	}
	return 1;
}
