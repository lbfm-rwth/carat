/* author: Oliver Heidbuechel */
/* last change: 23.11.2000 by author */

#include<typedef.h>
#include<datei.h>
#include<getput.h>
#include<matrix.h>
#include<longtools.h>
#include<tools.h>
#include"zass.h"
#include <bravais.h>
#include <presentation.h>
#include <graph.h>



/* --------------------------------------------------------------- */
/* constructs a matrix and writes M no-times on the diagonal       */
/* --------------------------------------------------------------- */
static matrix_TYP *to_diag(matrix_TYP *M,
                           int dim,
                           int no)
{
   int i, j, k;

   matrix_TYP *D;


   D = init_mat(M->rows, M->cols * no, "");

   for (i = 0; i < no; i++){
      for (j = 0; j < dim; j++){
         for (k = 0; k < dim; k++){
            D->array.SZ[i * dim + j][i * dim + k] = M->array.SZ[i * dim + j][k];
         }
      }
   }

   return(D);
}


/* --------------------------------------------------------------- */
/* transl_aff_normal                                               */
/* calculate a generating set of the linear part of                */
/* N_(Aff(R^n))(R(G,0)) with tranlations in Q^n/Z^n                */
/* returns the set {((G->gen[i] - I) * Translation_j)_i | j = ...} */
/* --------------------------------------------------------------- */
/* erzeuger: generators of the point group G                       */
/* erzanz  : number of matrices in "erzeuger"                      */
/* anzahl  : the number of the translations will be saved here     */
/* --------------------------------------------------------------- */
matrix_TYP **transl_aff_normal(matrix_TYP **erzeuger, 
                               int erzanz,
                               int *anzahl)
{
   matrix_TYP  *B, *B_diag,
              *translation,
              **tmp,
              **transl;
   
   int  k, j, i, dim, zahl;


   dim = erzeuger[0]->rows;

   /* create the matrix (g_1 - I, g_2 - I, ...)^tr */
   B = calc_B(erzeuger,erzanz);

   /* solve (g_j-I)s = 0 (Z^n) for all j */
   tmp = cong_solve(B);

   /* create matrices */
   for (anzahl[0] = 0;
        anzahl[0] < tmp[3]->cols && tmp[3]->array.SZ[anzahl[0]][anzahl[0]] != 0;
        anzahl[0]++);
    
   if (anzahl[0] != 0){
      /* create matrix diag (g_1 - I, g_2 - I, ...) */
      B_diag = to_diag(B, dim, erzanz);
      transl = (matrix_TYP **)calloc(anzahl[0], sizeof(matrix_TYP *));
      for (j = 0; j < anzahl[0]; j++){
         translation = init_mat(dim * erzanz, 1,"");
         for (i = 0; i < dim; i++){
            zahl = tmp[2]->array.SZ[i][j];
            for (k = 0; k < erzanz; k++){
               translation->array.SZ[i + k * dim][0] = zahl;
            }
         }
         translation->kgv = tmp[3]->array.SZ[j][j];
         translation->flags.Integral = FALSE;
         Check_mat(translation);
         transl[j] = mat_mul(B_diag, translation);
         Check_mat(transl[j]);
         if (transl[j]->kgv != 1){  /* paranoia test */
            fprintf(stderr, "ERROR in transl_aff_normal\n");
            exit(7);
         }
         free_mat(translation);
      }
      free_mat(B_diag);
   }
   else{
      transl = NULL;
   }   

   /* cleaning up */    
   for (j = 0; j < 4; j++){
      if (tmp[j] != NULL)
         free_mat(tmp[j]);
   }
   free(tmp); 
   free_mat(B);

   return(transl);
}



/* ---------------------------------------------------------------- */
/* cong_solve_part                                                  */
/* Calculates one solution of A*s = c (Z^m)                         */
/* If there is no solution, the function returns NULL.              */
/* ---------------------------------------------------------------- */
/* A and c are matrices with A in Z^(nxm) and c in Q^(nx1)          */
/* ---------------------------------------------------------------- */
static matrix_TYP *cong_solve_part(matrix_TYP *A,
                                   matrix_TYP *c)
{
  matrix_TYP *R,
             *L,
             *D,
             *r,
             *hilfsloesung,
             *loesung;

  int i, last;


  /* trivial case */
  if (equal_zero(A) == 1 && equal_zero(c) == 1){
     loesung = init_mat(A->cols, 1, "");
     return(loesung);
  }


  L = init_mat(A->rows,A->rows,"1");
  R = init_mat(A->cols,A->cols,"1");
  hilfsloesung = init_mat(A->cols,1,"");

  /* find matrices D,L,R with LAR = D and D is the elemantary divisor matrix
     and L in Gl_n(Z) and R in Gl_m(Z) */
  D = long_elt_mat(L, A, R);
  for (last = 0; last<D->cols && D->array.SZ[last][last] != 0; last++);
  r =  mat_mul(L,c);

  /* change the entries in r such that 0 <= r->array.SZ[i][0] < r->kgv
     (example: -5/3 will be changed to 1/3) */
  for (i=0; i<r->rows; i++){
     r->array.SZ[i][0] %= r->kgv;
     if (r->array.SZ[i][0] < 0)
        r->array.SZ[i][0] += r->kgv;
  }

  /* Find one solution of D*t = r (Z^m) */
  hilfsloesung->kgv = r->kgv * D->array.SZ[last-1][last-1];
  for (i=0; i<last; i++){
     if (r->array.SZ[i][0] == 0)
        hilfsloesung->array.SZ[i][0] = 0;
     else{
        hilfsloesung->array.SZ[i][0] = r->array.SZ[i][0]
               * D->array.SZ[last-1][last-1] / D->array.SZ[i][i];
     }
  }
  for (i=last; i<D->cols; i++){
     if (r->array.SZ[i][0] != 0){
        loesung = NULL;
     }
  }

  /* calculate one solution of A*s = c (Z^m) */
  loesung = mat_mul(R,hilfsloesung);

  /* cleaning up */
  free_mat(L);
  free_mat(R);
  free_mat(D);
  free_mat(r);
  free_mat(hilfsloesung);

  Check_mat(loesung);
  return(loesung);
}



/* --------------------------------------------------------------------- */
/* Let R be a spacegroup, P = P(R), coz = cocycle of R                   */
/* returns one element in the affine normalizer of R with linear part    */
/* given by 'lin' (has to be possible)                                   */
/* --------------------------------------------------------------------- */
matrix_TYP *to_aff_normal_element(matrix_TYP *lin,
                                  matrix_TYP *coz,
                                  int flag,
                                  bravais_TYP *P,
                                  bravais_TYP *R)
{
   matrix_TYP *N, *cozl, *DIAG, *cozr, *ccc, *B, *lin_inv, *part;

   bravais_TYP *conj;

   int d, j, k;

   rational eins, minuseins;



   /* prepare */
   eins.z = eins.n = minuseins.n = 1;
   minuseins.z = -1;
   d = lin->rows;
   N = copy_mat(lin);
   real_mat(N, d + 1, d + 1);
   N->array.SZ[d][d] = 1;

   if (flag != 1){
      conj = init_bravais(P->dim);
      conj->gen_no = P->gen_no;
      conj->gen = (matrix_TYP **)calloc(P->gen_no, sizeof(matrix_TYP *));
      lin_inv = mat_inv(lin);
      for (j = 0; j < P->gen_no; j++){
         conj->gen[j] = mat_kon(lin, P->gen[j], lin_inv);
      }
      free_mat(lin_inv);
      cozl = sg(R, conj);

      DIAG = matrix_on_diagonal(lin, P->gen_no);
      cozr = mat_mul(DIAG,coz);
      free_mat(DIAG);
      ccc = mat_add(cozl,cozr,minuseins,eins);
      free_mat(cozl);
      free_mat(cozr);
      Check_mat(ccc);
      B =  calc_B(conj->gen, conj->gen_no);
      part = cong_solve_part(B,ccc);
      if (part == NULL){
         printf("ERROR in aff_normalisator! cong_solve_part returns NULL!\n");
         exit(23);
      }
      free_mat(ccc);
      free_mat(B);
      free_bravais(conj);
   }
   else{
      /* one solution is zero if the spacegroup splits */
      part = init_mat(d, 1, "");
   }

   /* construct the matrix */
   if (part->kgv != 1){
      for (j = 0; j < d; j++)
         for (k = 0; k < d; k++)
            N->array.SZ[j][k] *= part->kgv;
      N->array.SZ[d][d] = part->kgv;
      N->kgv = part->kgv;
   }
   for (j = 0; j < d; j++)
      N->array.SZ[j][d] = part->array.SZ[j][0];

   Check_mat(N);
   free_mat(part);

   return(N);
}


