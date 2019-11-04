#include"typedef.h"
#include"matrix.h"
#include"symm.h"
#include "tools.h"
#include"bravais.h"

/****************************************************************************\
@  -----------------------------------------------------------------------
@  FILE: vor_neighbour.c:
@  -----------------------------------------------------------------------
@
\****************************************************************************/

/****************************************************************************\
@
@ matrix_TYP *voronoi_neighbour(A, X, Amin, lc, rc)
@ matrix_TYP *A, *X;
@ int Amin, *lc, *rc;
@
@ The arguments of 'voronoi_neigbour' are
@ A:        a positiv definite matrix
@ Amin:     The Minimum of A, min(A) :=  min{yAy^{tr} | y in Z^n, y not 0}
@ X:        a symmetric matrix, which is not positiv semidefinite 
@ lc, rc:   These pointer to an integer return the values of the
@           coefficients N = lc * A + rc *X for the result N of the
@           function
@
@  For A positive definite let M(A) denote the set of all y in Z^n
@  with yAy^{tr} = min(A) and
@  Z(A, X) the set of all y in M(A) with yXy^{tr} = 0.
@  The function voronoi_neigbbour calculates a Matrix N = lc *A + rc * X
@  with positive integral coefficients lc and rc such that
@  Z(A, X) is a proper subset of M(N).
@  In particular  min(N) = lc * min(A).
@-----------------------------------------------------------------------
@ 
\****************************************************************************/
matrix_TYP *
voronoi_neighbour (matrix_TYP *A, matrix_TYP *X, int Amin, int *lc, int *rc)
{
   int i,j,n;
   matrix_TYP *SV, *N;
   int lo,ro, lu, ru;
   int lm, rm;
   int anz, found;
   int merk, Xy, Ay;


   j = definite_test(X);
   if(j>= 0)
      return(NULL);
   n = A->rows;
   N = init_mat(n,n,"");
   N->flags.Integral = N->flags.Symmetric = TRUE;
   N->flags.Diagonal = FALSE;
  /************************************************************************\
  | First calculate Form N = lm * A + rm * X, with N positiv definite
  | and min(N) < lm * min(A)
  \************************************************************************/
  lu = 1; ru = 0;
  lo = 0; ro = 1;
  found = FALSE;
  while(found == FALSE)
  {
     lm = lu + lo;
     rm = ru + ro;
     for(i=0;i<n;i++)
       for(j=0;j<=i;j++)
         N->array.SZ[i][j]  = lm * A->array.SZ[i][j] + rm * X->array.SZ[i][j];
     for(i=0;i<n;i++)
      for(j=0;j<i;j++)
        N->array.SZ[j][i] = N->array.SZ[i][j]; 
     j = definite_test(N);
     if(j == 2)
     {
       SV = short_vectors(N, (lm * Amin -1), 0, 1, 0, &anz);
       if(anz > 0)
          found = TRUE;
       else
          { lu = lm; ru = rm; free_mat(SV);}
     }
     else
     { lo = lm; ro = rm;}
  }

  /************************************************************************\
  | The 1st row of SV containes an vector y with 
  |     lm * min(A) > yNy^{tr} = lm * yAy^{tr} + rm * yXy^{tr}
  |  Now yAy^{tr} > min(A) and yXy^{tr} < 0.
  |  Calculate the rational number p/q = (Amin - yAy^{tr}) / (yXy^{tr}) >= 0.
  |  Then y (q * A + p * X) y^{tr} = q * Amin
  |  Replace N by q * A + p X.
  |  If min(N) = q * min(A) we are done, otherwise caluclate vector y'
  |  with y'Ny'^{tr} < q * min(A) and repeat.
  \************************************************************************/
  found = FALSE;
  while(found == FALSE)
  {
    /************************************************************************\
    | calculate Ay = y A y^{tr}
    \************************************************************************/
     Ay = 0;
     for(i=0;i<n;i++)
     {
       merk = 0;
       for(j=0;j<n;j++)
         merk += A->array.SZ[i][j] * SV->array.SZ[0][j];
        Ay += SV->array.SZ[0][i] * merk;
     }
    /************************************************************************\
    | calculate Xy = y X y^{tr}
    \************************************************************************/
     Xy = 0;
     for(i=0;i<n;i++)
     {
       merk = 0;
       for(j=0;j<n;j++)
         merk += X->array.SZ[i][j] * SV->array.SZ[0][j];
        Xy += SV->array.SZ[0][i] * merk;
     }
    /************************************************************************\
    | calculate the new N
    \************************************************************************/
     Ay -= Amin;
     Xy = -Xy;
     merk = GGT(Ay, Xy);
     if(merk != 1)
     {
       Ay /= merk;
       Xy /= merk;
     }
     for(i=0;i<n;i++)
      for(j=0;j<=i;j++)
      {
        N->array.SZ[i][j] =  Xy * A->array.SZ[i][j] + Ay * X->array.SZ[i][j];
        N->array.SZ[j][i] = N->array.SZ[i][j];
      }
     free_mat(SV);
     SV = short_vectors(N, (Xy * Amin - 1), 0, 1, 0, &anz);
     if(anz == 0)
     { free_mat(SV); found = TRUE; *lc = Xy; *rc = Ay;}
  }

  /* canceled these lines on 16/2/97 because they don't fit
     to specification and cause problems: tilman
  ********************************************************************\
  | calculate the gcd of the entries of N and divide N by it
  \********************************************************************
  g = N->array.SZ[0][0];
  for(i=1;i<n;i++)
   for(j=0;j<=i;j++)
     g = GGT(g, N->array.SZ[i][j]);
  if(g < 0)
   g = -g;
  if(g != 1)
  {
    for(i=0;i<n;i++)
     for(j=0;j<n;j++)
      N->array.SZ[i][j] /= g;
  }

  printf("in vor_nei g %d\n",g);
  */

  return(N);
}
