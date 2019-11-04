#include "typedef.h"
#include "matrix.h"
#include "symm.h"
#include "bravais.h"
#include "voronoi.h"
#include "longtools.h"
#include "tools.h"
#include "getput.h"


/*********************************************************************\
@
@---------------------------------------------------------------------
@  FILE: first_perfect.c
@---------------------------------------------------------------------
@
\*********************************************************************/

/*********************************************************************\
@
@ matrix_TYP *first_perfect(A, G, Ftr, trbifo, min)
@ matrix_TYP *A, **Ftr, *trbifo;
@ bravais_TYP *G;
@ int *min;
@
@ 'first_perfect' calculates a G-perfect Form.
@ The arguments of 'first_perfect' are
@ A:		a G-invariant, positive definite Form
@ G:		the group for which a perfect form will be calculated
@ Ftr:		a Z-basis for the integral forms in the formspace of G^{tr}
@ trbifo:	The matrix of the bilinear form
@		F(G)xF(G^{tr})--> R: (A,B)---> trace(AB)
@		with respect to the Z-bases given by G->form and Ftr
@ min:		The pointer min returns the minimum of the calculated
@               perfect form.
@---------------------------------------------------------------------
@
\*********************************************************************/
matrix_TYP *
first_perfect (matrix_TYP *A, bravais_TYP *G, matrix_TYP **Ftr, matrix_TYP *trbifo, int *min)
{
  int i,j,k, n, lc, rc, g;
  matrix_TYP **V, *Vmat, *P, *M, *W, *W1, *M1;
  int Vanz;
  int rang;
  int Pmin;


  n = G->form_no;
  P = copy_mat(A);
  /*****************************************************************\
  | Check whether P is G-perfect or not
  \*****************************************************************/
  V = voronoi_vertices(P, G, &Vanz, &Pmin, &i);
  Vmat = init_mat(Vanz, G->form_no, "");
  for(i=0;i<Vanz;i++)
  {
    form_to_vec(Vmat->array.SZ[i], V[i], Ftr, n, &k);
    free_mat(V[i]);
  }
  free(V);
  rang = long_row_basis(Vmat,FALSE);
  if(rang == n)
  {
     free_mat(Vmat);
     Vmat = NULL;
     /* inserted 11/2/97 tilman */
     Vmat = shortest(P,min);
     free_mat(Vmat);

     return(P);
  }
  W = init_mat(n, n, "");
  M = init_mat(A->rows, A->cols, "");
  while(rang != n)
  {
     /* printf("The Voronoi domain of the following form has dimension %d\n", rang);
     put_mat(P, NULL, "next approximation", 2); */
     /****************************************************************\
     | Search for a form M perpendicular to the Voronoi-vertices
     |  (with respect to trbifo)
     | if M is positive semidefinite take -M
     \****************************************************************/
     for(i=0;i<rang;i++)
       for(j=0;j<n;j++)
       {
         W->array.SZ[i][j] = 0;
          for(k=0;k<n;k++)
           W->array.SZ[i][j] += Vmat->array.SZ[i][k] * trbifo->array.SZ[j][k];
       }
     free_mat(Vmat);
     Vmat = NULL;
     W1 = long_kernel_mat(W);
     for(i=0;i<A->rows;i++)
       for(j=0;j<=i;j++)
       {
         M->array.SZ[i][j] = 0;
         for(k=0;k<n;k++)
           M->array.SZ[i][j] += W1->array.SZ[k][W1->cols-1] * G->form[k]->array.SZ[i][j];
         M->array.SZ[j][i] = M->array.SZ[i][j];
       }
     free_mat(W1);
     k = definite_test(M);
     if(k >= 0)
     {
        for(i=0;i<M->rows;i++)
          for(j=0;j<M->cols;j++)
            M->array.SZ[i][j] = -M->array.SZ[i][j];
     }
     /*****************************************************************\
     | calculate next form with more shortest vectors
     | and divide the matrix by the gcd of the entries
     \*****************************************************************/
     M1 = voronoi_neighbour(P, M, Pmin, &lc, &rc);
     free_mat(P);
     P = M1;
     M1 = NULL;
     Pmin = Pmin * lc;
     g = Pmin;
     for(i=0;i<P->rows && g != 1; i++)
       for(j=0;j<=i && g != 1;j++)
       {
         if(P->array.SZ[i][j] != 0)
         {  k = GGT(g, P->array.SZ[i][j]); g = k;}
       }
     if(g != 1 && g!= -1)
     {
        Pmin /= g;
        for(i=0;i<P->cols;i++)
         for(j=0;j<=i;j++)
         {
             P->array.SZ[i][j] /= g;
             P->array.SZ[j][i] = P->array.SZ[i][j];
         }
     }

     /*****************************************************************\
     | if rang < n-1, check whether P is G-perfect or not
     | if rang = n-1, the new Form must be perfect
     \*****************************************************************/
     if(rang < n-1)
     {
       V = voronoi_vertices(P, G, &Vanz, &k, &i);
       if(k != Pmin)
       {
         put_mat(P,NULL,"P",2);
         printf("Minimum of P %d, Pmin %d\n",k,Pmin);
         printf("Fehler in first_perfect, Minimum falsch\n");
         exit(3);
       }
       Vmat = init_mat(Vanz, G->form_no, "");
       for(i=0;i<Vanz;i++)
       {
         form_to_vec(Vmat->array.SZ[i], V[i], Ftr, n, &k);
         free_mat(V[i]);
       }
       free(V);
       rang = long_row_basis(Vmat,FALSE);
     }
     else{
       rang = n;
     }
  }
  free_mat(W);
  free_mat(M);
  *min = Pmin;

  if (Vmat != NULL){
     free_mat(Vmat);
  }

  return(P);
}
