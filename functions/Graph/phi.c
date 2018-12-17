/* last change: 02.11.00 by Oliver Heidbuechel */



#include <typedef.h>
#include <matrix.h>
#include <bravais.h>
#include <base.h>
#include <graph.h>
#include <zass.h>
#include <datei.h>
#include <longtools.h>
#include <tools.h>


/* ----------------------------------------------------------------------------- */
/* this is a function which make the same as the program same_generators         */
/* the code is taken out of the main-function of same_generators                 */
/* ----------------------------------------------------------------------------- */
matrix_TYP *sg(bravais_TYP *R,
               bravais_TYP *P)
{
   matrix_TYP **RG, *coz;

   int **words, denominator,
       j, k;



   RG = (matrix_TYP **)calloc(R->gen_no, sizeof(matrix_TYP *));
   words = (int **)calloc(P->gen_no, sizeof(int *));
   denominator = 1;
   for (j = 0; j < R->gen_no; j++){
      RG[j] = copy_mat(R->gen[j]);
      RG[j]->cols--;
      RG[j]->rows--;
      Check_mat(RG[j]);
      if (!RG[j]->flags.Integral){
         fprintf(stderr,"The point group has to be integral\n");
         exit(3);
      }
      rat2kgv(R->gen[j]);
      denominator *= (R->gen[j]->kgv / GGT(R->gen[j]->kgv, denominator));
   }

   /* stick the rigth INTEGRAL cozycle at the end of the RG[j] */
   for (j = 0; j < R->gen_no; j++){
      RG[j]->cols++;
      RG[j]->rows++;
      for (k = 0; k < RG[j]->rows-1; k++){
         RG[j]->array.SZ[k][R->dim-1] = (denominator / R->gen[j]->kgv) *
                   R->gen[j]->array.SZ[k][R->dim-1];
         RG[j]->array.SZ[R->dim-1][R->dim-1] = 1;
         Check_mat(RG[j]);
      }
   }

   /* get the cozycle on the right generators */
   coz = reget_gen(RG, R->gen_no, P, words, TRUE);

   /* the cozykle has to become the right denominator */
   coz->kgv = denominator;
   Check_mat(coz);

   for (j = 0; j < R->gen_no; j++){
      free_mat(RG[j]);
   }
   for (j = 0; j < P->gen_no; j++){
      free(words[j]);
   }
   free(words);
   free(RG);

   return(coz);
}



/* ----------------------------------------------------------------------------- */
/* GL = <gl_1, ..., gl_m> = <s1, ..., s_m> = standard                            */
/* X: informations about the cohomology group of standard                        */
/* calculate the cocycles for g1_1, ..., gl_m                                    */
/* ----------------------------------------------------------------------------- */
matrix_TYP *H1_of_standard_to_GL(bravais_TYP *GL,
                                 bravais_TYP *standard,
                                 matrix_TYP **X)
{
   int i, j, X_first;

   matrix_TYP *coz, *help;

   bravais_TYP *RG;



   for (X_first = 0; X_first < X[1]->cols && X[1]->array.SZ[X_first][X_first] == 1; X_first++);
   coz = init_mat(X[0]->rows, X[0]->cols, "");

   for (i = 0; i < X[0]->cols; i++){
      help = init_mat(X[0]->rows, 1, "");
      for (j = 0; j < X[0]->rows; j++){
         help->array.SZ[j][0] = X[0]->array.SZ[j][i];
      }
      help->kgv = X[1]->array.SZ[i + X_first][i + X_first];
      help->flags.Integral = FALSE;
      RG = extract_r(standard, help);
      free_mat(help);
      help = sg(RG, GL);
      if (help->kgv != X[1]->array.SZ[i + X_first][i + X_first]){
         fprintf(stderr, "ERROR in H1_of_standard_to_GL!\n");
         exit(17);
      }
      for (j = 0; j < X[0]->rows; j++){
         coz->array.SZ[j][i] = help->array.SZ[j][0];
      }
      free_mat(help);
      free_bravais(RG);
   }

   return(coz);
}



/* ----------------------------------------------------------------------------- */
/* calculate the kernel and H^1/Ker for phi                                      */
/* ----------------------------------------------------------------------------- */
static void kernel_fkt(matrix_TYP *phi,
                       matrix_TYP *A,
                       int A_first,
                       matrix_TYP *B,
                       matrix_TYP **kernel,
                       H1_mod_ker_TYP *H1_mod_ker)
{
   int i, j, rows, cols, B_first, counter = 0, flagge, zahl,
       D_first, D_last, D_diff;

   matrix_TYP *ker, *L, *R, *D, *Li;


   for (B_first = 0; B_first < B->cols && B->array.SZ[B_first][B_first] == 1; B_first++);

   /* expand phi */
   rows = phi->rows;
   cols = phi->cols;
   real_mat(phi, rows, cols + rows);
   for (i = 0; i < rows; i++){
      phi->array.SZ[i][i + cols] = B->array.SZ[i + B_first][i + B_first];
   }

   /* calculate kernel */
   ker = long_kernel_mat(phi);
   kernel[0] = init_mat(cols, ker->cols, "");
   for (i = 0; i < ker->cols; i++){
      flagge = 0;
      for (j = 0; j < cols; j++){
         zahl = ker->array.SZ[j][i] % A->array.SZ[j + A_first][j + A_first];
         if (zahl < 0)
            zahl += A->array.SZ[j + A_first][j + A_first];
         if (zahl != 0)
            flagge = 1;
         kernel[0]->array.SZ[j][counter] = zahl;
      }
      if (flagge == 1){
         counter++;
      }
   }
   free_mat(ker);

   /* calculate A/ker(phi) */
   if (counter == 0){
      /* trivial part: ker(phi) = 0 */
      H1_mod_ker[0].D = init_mat(cols, cols, "");
      H1_mod_ker[0].M = (matrix_TYP **)calloc(cols, sizeof(matrix_TYP *));
      H1_mod_ker[0].i = init_mat(cols, cols, "1");
      for (i = 0; i < cols; i++){
         H1_mod_ker[0].D->array.SZ[i][i] = A->array.SZ[i + A_first][i + A_first];
         H1_mod_ker[0].M[i] = init_mat(cols, 1, "");
         H1_mod_ker[0].M[i]->array.SZ[i][0] = 1;
      }
      H1_mod_ker[0].erz_no = cols;
      H1_mod_ker[0].D_first = 0;
   }
   else{
      real_mat(kernel[0], cols, counter + cols);
      for (i = 0; i < cols; i++){
         kernel[0]->array.SZ[i][i + counter] = A->array.SZ[i + A_first][i + A_first];
      }
      L = init_mat(cols, cols, "1");
      R = init_mat(counter + cols, counter + cols, "1");
      D = long_elt_mat(L, kernel[0], R);
      Li = mat_inv(L);
      /*
      standard_form(L, A, A_first);
      Li = graph_mat_inv(L, A, A_first);
      */
      for (D_first = 0; D_first < D->rows && D->array.SZ[D_first][D_first] == 1; D_first++);
      for (D_last = D_first; D_last < D->rows && D->array.SZ[D_last][D_last] != 0; D_last++);
      D_diff = D_last - D_first;
      H1_mod_ker[0].M = (matrix_TYP **)calloc(D_diff, sizeof(matrix_TYP *));
      for (i = 0; i < D_diff; i++){
         H1_mod_ker[0].M[i] = init_mat(cols, 1, "");
         for (j = 0; j < cols; j++){
            H1_mod_ker[0].M[i]->array.SZ[j][0] = Li->array.SZ[j][i + D_first] %
                                              A->array.SZ[j + A_first][j + A_first];
            if (H1_mod_ker[0].M[i]->array.SZ[j][0] < 0)
               H1_mod_ker[0].M[i]->array.SZ[j][0] += A->array.SZ[j + A_first][j + A_first];
         }
      }
      H1_mod_ker[0].erz_no = D_diff;
      H1_mod_ker[0].D = D;
      H1_mod_ker[0].D_first = D_first;

      H1_mod_ker[0].i = L;
      free_mat(Li);
      free_mat(R);
      D = NULL;
      Li = NULL;
   }

   /* clean */
   real_mat(phi, rows, cols);
   real_mat(kernel[0], cols, counter);
}



/* -------------------------------------------------------------------------------- */
/* calculate phi: H^1(G, Q^n/L) -> H^1(G,Q^n/Z^n), i.e.                             */
/* calculate phi: H^1(given by Xi) -> H^1(given by Xj)                              */
/* -------------------------------------------------------------------------------- */
void calculate_phi(matrix_TYP *diag,
                   matrix_TYP *coz,
                   matrix_TYP **Xi,
                   matrix_TYP **Xj,
                   matrix_TYP *GLS,
                   matrix_TYP **phi,
                   matrix_TYP **kernel,
                   matrix_TYP ***image,
                   int *image_gen_no,
                   H1_mod_ker_TYP *H1_mod_ker)
{
   matrix_TYP *H1_L,
              *tmp,
*test;

   int i, j,
       dim_i, dim_j, first;



   for (first = 0; first < Xj[1]->cols && Xj[1]->array.SZ[first][first] == 1; first++);

   /* calculate cohomology group for the lattice */
   if (coz->cols > 0)
      H1_L = mat_mul(diag, coz);
   else
      H1_L = init_mat(diag->rows, 0, "");

   dim_j = H1_L->cols;
   dim_i = Xi[0]->cols;
   phi[0] = init_mat(dim_i, dim_j, "");

   /* calculate phi, kernel and image */
   if (dim_i > 0 && dim_j > 0){
      image[0] = (matrix_TYP **)calloc(dim_j, sizeof(matrix_TYP *));
      /* some functions need the generating set for the image in the following form */
      for (i = 0; i < dim_j; i++){
         tmp = init_mat(H1_L->rows, 1, "");
         for (j = 0; j < tmp->rows; j++)
            tmp->array.SZ[j][0] = H1_L->array.SZ[j][i];
         tmp->kgv = Xj[1]->array.SZ[i + first][i + first];
         tmp->flags.Integral = 0;
         Check_mat(tmp);
         image[0][i] = standard_rep(tmp, GLS, Xi[1]);
         free_mat(tmp);
         for (j = 0; j < dim_i; j++){
            phi[0]->array.SZ[j][i] = image[0][i]->array.SZ[j][0];
         }
      }
      kernel_fkt(phi[0], Xj[1], first, Xi[1], kernel, H1_mod_ker);
      image_gen_no[0] = dim_j;
      H1_mod_ker[0].flag = 0;
   }
   else{
      /* trivial part */
      if (dim_i == 0){
         kernel[0] = init_mat(dim_j, dim_j, "1");
      }
      else{
         kernel[0] = init_mat(dim_j, 0, "");
      }
      image[0] = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
      image_gen_no[0] = 0;
      if (dim_i == 0 && dim_j == 0)
         H1_mod_ker[0].flag = 3;
      else{
         if (dim_i == 0)
            H1_mod_ker[0].flag = 2;
         else
            H1_mod_ker[0].flag = 1;
      }
   }

   /* clean */
   free_mat(H1_L);

}




