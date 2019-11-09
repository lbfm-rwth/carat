#include "typedef.h"
#include "matrix.h"
#include "voronoi.h"
#include "ZZ_P.h"

static matrix_TYP *Mod, *Trf;
static rational c;

/*{{{}}} */
/*{{{  lll_init, static */
static void 
lll_init (matrix_TYP *Mat, int lll_bnd)

{
  int **Z;
  int i;


  c.z = lll_bnd;
  c.n = lll_bnd + 1;

  Trf = mat_red (Mat);
  dec_mat (Mat, Trf);

  if (LLLREDUCED)
    {
      Mod = init_mat (Mat->rows, Mat->rows, "srl");
      Z = Mod->array.SZ;
      for (i = 0; i < Mat->rows; i++)
	{
	  Z[i][i] = 1;
	}
    }
}
/*}}}  */
/*{{{  upd_mod, static */

static void 
upd_mod (matrix_TYP *Mat, int pos)

{
  int **Z, **N, **A;
  rational sum, help;
  int i, j;

  sum = help = Zero;

  A = Mat->array.SZ;
  Z = Mod->array.SZ;
  N = Mod->array.N;
  if (pos == 1)
    Z[0][0] = A[0][0];
  for (i = 0; i <= pos; i++)
    {
      sum.z = 0;
      sum.n = 1;
      for (j = 0; j < i; j++)
	{
	  help.z = Z[j][j] * Z[pos][j] * Z[i][j];
	  help.n = N[j][j] * N[pos][j] * N[i][j];
	  sum.z *= help.n;
	  help.z *= sum.n;
	  sum.z += help.z;
	  sum.n *= help.n;
	  Normal (&sum);
	}
      if ((i != pos) && (Z[i][i] != 0))
	{
	  N[pos][i] = sum.n * Z[i][i];
	  sum.n = sum.n * A[pos][i] - sum.z;
	  Z[pos][i] = sum.n * N[i][i];
	}
      else
	{
	  N[pos][i] = sum.n;
	  sum.n *= A[pos][i];
	  Z[pos][i] = sum.n - sum.z;
	}
      Normal2 (&Z[pos][i], &N[pos][i]);
    }
}

/*}}}  */
/*{{{  reduce, static */
static boolean 
reduce (matrix_TYP *Mat, int k, int l)
{
  rational d;
  int **Z, **N, **A, **T, *T_help;
  int i;
  int f;

  d = Zero;

  Z = Mod->array.SZ;
  N = Mod->array.N;
  A = Mat->array.SZ;
  T = Trf->array.SZ;
  d.z = abs (Z[k][l]);
  d.n = N[k][l];
  if ((d.z * 11) > d.n)
    {
      f = d.z / d.n;
      d.z %= d.n;
      if (d.z > (d.n / 2))
	f++;
      if (f != 0)
	{
	  if (Z[k][l] < 0)
	    f *= -1;
	  for (i = 0; i < l; i++)
	    {
	      Z[k][i] =
		Z[k][i] * N[l][i] - f * N[k][i] * Z[l][i];
	      N[k][i] *= N[l][i];
	      Normal2 (&Z[k][i], &N[k][i]);
	    }
	  Z[k][i] -= f * N[k][i];
	  A[k][k] -= (2 * A[k][l] - f * A[l][l]) * f;
	  for (i = 0; i <= l; i++)
	    A[k][i] -= f * A[l][i];
	  for (; i < k; i++)
	    A[k][i] -= f * A[i][l];
	  for (i++; i < Mat->rows; i++)
	    A[i][k] -= f * A[i][l];
	  for (i = 0; i < Trf->cols; i++)
	    T[k][i] -= f * T[l][i];
	}
      if (A[k][k] == 0)
	{
	  kill_row (Mat, k);
	  kill_col (Mat, k);
	  T_help = T[k];
	  /* Mat->cols = --Mat->rows; */
	  for (i = k; i < Mat->rows; i++)
	    {
	      T[i] = T[i + 1];
	      /* A[i] = A[i+1];
	         for ( j = k; j <= i; j++)
	         A[i][j]= A[i][j+1]; A[i][j+1] = 0; */
	    }
	  T[i] = T_help;
	  A = Mat->array.SZ;
	  return (FALSE);
	}
    }
  return (TRUE);
}
/*}}}  */
/*{{{  swap, static */

static void 
swap (matrix_TYP *Mat, int i)

{
  rational help;
  int **Z, **N, **A, **T, *t;
  int j;
  int k;

  help = Zero;

  A = Mat->array.SZ;
  Z = Mod->array.SZ;
  N = Mod->array.N;
  T = Trf->array.SZ;
  t = T[i + 1];
  T[i + 1] = T[i];
  T[i] = t;
  for (j = 0; j < i; j++)
    {
      k = A[i + 1][j];
      A[i + 1][j] = A[i][j];
      A[i][j] = k;
      Z[i][j] = Z[i + 1][j];
      Z[i + 1][j] = 0;
      N[i][j] = N[i + 1][j];
      N[i + 1][j] = 1;
    }
  k = A[i][i];
  A[i][i] = A[i + 1][i + 1];
  A[i + 1][i + 1] = k;
  for (j += 2; j < Mat->rows; j++)
    {
      k = A[j][i + 1];
      A[j][i + 1] = A[j][i];
      A[j][i] = k;
    }
  help.z = Z[i + 1][i] * Z[i + 1][i] * Z[i][i] * N[i + 1][i + 1];
  help.n = N[i + 1][i] * N[i + 1][i] * N[i][i];
  Normal (&help);
  N[i][i] = help.n * N[i + 1][i + 1];
  Z[i][i] = help.n * Z[i + 1][i + 1] + help.z;
  Normal2 (&Z[i][i], &N[i][i]);
}

/*}}}  */
/*{{{  condition, static */
static boolean 
condition (int i)
{
  rational help;
  int **Z, **N;
  boolean cond;

  help = Zero;
  Z = Mod->array.SZ;
  N = Mod->array.N;

  help.z = Z[i][i - 1] * Z[i][i - 1];
  help.n = N[i][i - 1] * N[i][i - 1];
  help.z = -help.z * c.n + c.z * help.n;
  help.n *= c.n;
  help.z *= Z[i - 1][i - 1];
  help.n *= N[i - 1][i - 1];
  help.n *= Z[i][i];
  help.z *= N[i][i];
  cond = help.n < help.z;
  return (cond);
}

/*}}}  */
/*{{{  lll, exported */
matrix_TYP *
ZZ_lll (matrix_TYP *Mat, int lll_bnd)

{
  int i = 1, j;

  lll_init (Mat, lll_bnd);
  if (LLLREDUCED)
    {
      if (Mat->cols == 1)
	return (Trf);
      do
	{
	  upd_mod (Mat, i);
	  for (j = i - 1; j >= 0; j--)
	    {
	      if (!reduce (Mat, i, j))
		{
		  i--;
		  break;
		}

	      if ((j == i - 1) && condition (i))
		{
		  swap (Mat, --i);
		}
	    }
	  i++;
	} while (i < Mat->cols);
      free_mat (Mod);
      for (i = 0; i < Mat->rows; i++)
	for (j = i + 1; j < Mat->cols; j++)
	  Mat->array.SZ[i][j] = Mat->array.SZ[j][i];
    }
  Check_mat (Trf);
  Check_mat (Mat);
  return (Trf);
}
/*}}}  */
