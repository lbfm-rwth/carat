#include "typedef.h"
#include "matrix.h"
#include "bravais.h"
#include "longtools.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: lincomb.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *vec_to_form(v, F, Fanz)
@ int *v;
@ matrix_TYP **F;
@ int Fanz;
@
@ calculates the matrix v[1] F[1]  + v[2] F[2] +....+v[Fanz-1] F[Fanz-1}
@---------------------------------------------------------------------------
@
\**************************************************************************/

matrix_TYP *vec_to_form(v, F, Fanz)
int *v;
matrix_TYP **F;
int Fanz;
{
   int i,j,k, n,m;
   matrix_TYP *M;
   int sym;

   n = F[0]->rows;
   m = F[0]->cols;
   for(i=1;i<Fanz;i++)
   {
       if(F[i]->rows != n || F[i]->cols != m)
       {
          printf(" Different size of matrices in vec_to_form\n");
          exit(3);
       }
   }
   i=0; sym = TRUE;
   while(i<Fanz && F[i]->flags.Symmetric == TRUE)
     i++;
   if(i<Fanz)
      sym = FALSE;

  if(sym == FALSE)
  {
   M = init_mat(n,m, "");
   for(i=0;i<Fanz;i++)
   {
     if(v[i] != 0)
     {
       for(j=0;j<n;j++)
         for(k=0;k<m;k++)
           M->array.SZ[j][k] += v[i] * F[i]->array.SZ[j][k];
     }
   }
  }
  else
  {
   M = init_mat(n,m, "s");
   for(i=0;i<Fanz;i++)
   {
     if(v[i] != 0)
     {
       for(j=0;j<n;j++)
         for(k=0; k<=j; k++)
           M->array.SZ[j][k] += v[i] * F[i]->array.SZ[j][k];
     }
   }
   for(j=0;j<n;j++)
    for(k=0;k<j;k++)
      M->array.SZ[k][j] = M->array.SZ[j][k];
  }
  Check_mat(M);
  return(M);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ void form_to_vec(erg, A, F, Fanz, denominator)
@ int *erg;
@ matrix_TYP *A, **F;
@ int Fanz, *denominator;
@
@ calculates a vector v such that A = 1/d(v[0]F[0] +...+v[Fanz-1]F[Fanz-1]
@ and writes it onto the vector erg.
@ the space for erg must be allocated before using this function.
@ the integer d is returned via the pointer denominator.
@---------------------------------------------------------------------------
@
\**************************************************************************/
void form_to_vec(erg, A, F, Fanz, denominator)
int *erg;
matrix_TYP *A, **F;
int Fanz, *denominator;
{
  int i,j,k,l, n,m,r;
  matrix_TYP *X, *X1, *Trf;
  int sym, rang;

   n = F[0]->rows;
   m = F[0]->cols;
   for(i=1;i<Fanz;i++)
   {
       if(F[i]->rows != n || F[i]->cols != m)
       {
          printf(" Different size of matrices in form_to_vec\n");
          exit(3);
       }
   }
   i=0; sym = TRUE;
   while(i<Fanz && F[i]->flags.Symmetric == TRUE)
     i++;
   if(i<Fanz)
      sym = FALSE;

   if(sym == FALSE)
   {
     r = n*m;
     X = init_mat(r,Fanz+1, "");
     for(i=0;i<Fanz;i++)
     {
        l = 0;
        for(j=0;j<n;j++)
         for(k=0;k<n;k++)
         {  X->array.SZ[l][i] = F[i]->array.SZ[j][k]; l++;}
     }
     l = 0;
     for(j=0;j<n;j++)
      for(k=0;k<n;k++)
      {  X->array.SZ[l][Fanz] = A->array.SZ[j][k]; l++;}
   }
   else
   {
     r = (n*(n+1))/2;
     X = init_mat(r,Fanz+1, "");
     for(i=0;i<Fanz;i++)
     {
        l = 0;
        for(j=0;j<n;j++)
         for(k=0;k<=j;k++)
         {  X->array.SZ[l][i] = F[i]->array.SZ[j][k]; l++;}
     }
     l = 0;
     for(j=0;j<n;j++)
      for(k=0;k<=j;k++)
      {  X->array.SZ[l][Fanz] = A->array.SZ[j][k]; l++;}
   }
   X->rows = Fanz-1;
   do
   {
      X->rows++;

      /* tilman: changed 02/09/97 from (to avoid overflow)
      rang = long_row_gauss(X); */
      rang = long_row_basis(X,FALSE);

   }while(X->rows < r && rang != Fanz);

   if(X->rows == r && rang != Fanz)
   {
     printf("Matrices F in form_to_vec are not linear independent\n");
     exit(3);
   }
   X->rows = rang;
   X1 = tr_pose(X);
   X->rows = r;
   free_mat(X);
   Trf = init_mat(X1->rows, X1->rows, "");
   long_row_trf_gauss(X1, Trf);
   free_mat(X1);
   if(Trf->array.SZ[Fanz][Fanz] == 0)
   {
     printf("Matrices F in form_to_vec are not linear independent\n");
     exit(3);
   }
   if(Trf->array.SZ[Fanz][Fanz] < 0)
   {
     for(i=0;i<Fanz;i++)
       erg[i] = Trf->array.SZ[Fanz][i];
     *denominator = - Trf->array.SZ[Fanz][Fanz];
   }
   if(Trf->array.SZ[Fanz][Fanz] > 0)
   {
     for(i=0;i<Fanz;i++)
       erg[i] = - Trf->array.SZ[Fanz][i];
     *denominator = Trf->array.SZ[Fanz][Fanz];
   }
   free_mat(Trf);
}

/**************************************************************************\
@---------------------------------------------------------------------------
@ vertex_TYP *form_to_vertex(A, F, Fanz, denominator)
@ matrix_TYP *A, **F;
@
@ calculates a vector v such that A = 1/d(v[0]F[0] +...+v[Fanz-1]F[Fanz-1]
@ and returns it as V->v in a vertex_TYP V..
@ The integer d is returned via the pointer denominator.
@---------------------------------------------------------------------------
@
\**************************************************************************/
vertex_TYP *form_to_vertex(A, F, Fanz, denominator)
matrix_TYP *A, **F;
int Fanz, *denominator;
{
  vertex_TYP *v;
  extern vertex_TYP *init_vertex();
  v = init_vertex(Fanz, 0);
  form_to_vec(v->v, A, F, Fanz, denominator);
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ void form_to_vec_modular(erg, A, F, Fanz)
@ matrix_TYP *A, **F;
@ int *erg, Fanz;
@
@ the same as form_to_vec, but works only if the linear combination is
@ integral, i.e. the denominator is assumed to be 1.
@ the function calculated the result modulo big primes and fits it
@ together with the chinese remainder theorem.
@ This is faster as form_to_vec and avoids overflow.
@---------------------------------------------------------------------------
@
\**************************************************************************/
void form_to_vec_modular(erg, A, F, Fanz)
matrix_TYP *A, **F;
int *erg, Fanz;
{

  int i,j,k,l;
  int p1, p2, a1, a2, y1, y2;
  matrix_TYP **X1, **X2;
  matrix_TYP *G1, *G2;
  int X1anz, X2anz;
  int n, m;
  double q, res, q1, q2;
  double waste;
  int symmetric;

  p1 = 65521; p2 = 65519;
  n = A->cols;
  for(i=0;i<Fanz;i++)
  {
    if(F[i]->cols != n)
    {
     printf("error: matrices in 'form_to_vec_modular' have different number of cols\n");
     exit(3);
    }
  }
  n = A->rows;
  for(i=0;i<Fanz;i++)
  {
    if(F[i]->rows != n)
    {
     printf("error: matrices in 'form_to_vec_modular' have different number of rows\n");
     exit(3);
    }
  }

  symmetric = A->flags.Symmetric;
  for(i=0;i<Fanz;i++)
  {
     if(F[i]->flags.Symmetric == FALSE)
       symmetric = FALSE;
  }
  if(symmetric)
   n = (A->cols * (A->cols + 1))/2;
  else
   n = A->cols * A->rows;
  G1 = init_mat(n, Fanz, "");
  G2 = init_mat(n, 1, "");
  if(symmetric)
  {
    l = 0;
    for(i=0;i<A->cols;i++)
     for(j=0;j<=i;j++)
     {
       for(k=0;k<Fanz;k++)
         G1->array.SZ[l][k] = F[k]->array.SZ[i][j];
       G2->array.SZ[l][0] = A->array.SZ[i][j];
       l++;
     }
  }
  else
  {
    l = 0;
    for(i=0;i<A->rows;i++)
     for(j=0;j<A->cols;j++)
     {
       for(k=0;k<Fanz;k++)
         G1->array.SZ[l][k] = F[k]->array.SZ[i][j];
       G2->array.SZ[l][0] = A->array.SZ[i][j];
       l++;
     }
  }

  X1 = p_lse_solve(G1, G2, &X1anz, p1);
  X2 = p_lse_solve(G1, G2, &X2anz, p2);
  if(X1anz == 0 || X2anz == 0)
  {
     printf("error in 'form_to_vec_modular':\n");
     printf("The matrix A is not in the Z-span of the matrices F[0],..,F[%d],\n", (Fanz-1));
     exit(3);
  }
  if(X1anz != 1 || X2anz != 1)
  {
     printf("error in 'form_to_vec_modular':\n");
     printf("The matrix A is not in the Z-span of the matrices F[0],..,F[%d],\n", (Fanz-1));
     printf("or F[0],...,F[%d] are not independent over the field with", (Fanz-1));
     if(X1anz != 1)  printf("  %d  ", p1);
     if(X2anz != 1)  printf("  %d  ", p2);
     printf("elements\n");
     exit(3);
  }
  a1 = p_inv(p2, p1);
  a2 = p_inv(p1, p2);
  q = ((double) p1) * ((double) p2);
  q1 = ((double) p2) * ((double) a1);
  q2 = ((double) p1) * ((double) a2);

  for(i=0;i<Fanz;i++)
  {
    erg[i] = 0;
    y1 = X1[0]->array.SZ[i][0];
    y2 = X2[0]->array.SZ[i][0];
    if(y1 == y2)
      erg[i] = y1;
    else
    {
      res = ((double) y1) * q1 + ((double) y2) * q2;
      modf(res/q, &waste);
      res = res - waste * q;
      while( (res/q) > 0.5)
         res = res -q;
      while( (res/q) < -0.5)
         res = res + q;
      erg[i] = ((int) res);
    }
  }
  /*****************************************************************\
  | Ckeck if the modular result is also integrally correct
  \*****************************************************************/
  for(i=0; i<n;i++)
  {
     y1 = -(G2->array.SZ[i][0]);
     for(k=0;k<Fanz;k++)
       y1 += G1->array.SZ[i][k] * erg[k];
     if(y1 != 0)
     {
      printf("Sorry: Overflow in 'form_to_vec_modular'\n");
      exit(3);
     }
  }
  free_mat(X1[0]); free_mat(X2[0]);
  free(X1); free(X2);
  free_mat(G1); free_mat(G2);
}
