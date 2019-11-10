/* last change: 14.09.2000 by Oliver Heidbuechel */


#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <gmp.h>
#include <zass.h>
#include <longtools.h>
#include <orbit.h>
#include <bravais.h>
#include <graph.h>



/* -------------------------------------------------------------------- */
/* matrix_on_diagonal:							*/
/* Creates the matrix M = diag(mat, ... , mat)				*/
/* -------------------------------------------------------------------- */
matrix_TYP *matrix_on_diagonal(matrix_TYP *mat,
                               int anz)
{
  int i, j, k;

  matrix_TYP *M;


  M = init_mat(mat->rows * anz, mat->cols * anz, "");
  for (i = 0; i < anz; i++){
     for (j = 0; j < mat->rows; j++){
        for (k = 0; k < mat->cols; k++){
           M->array.SZ[i*mat->rows + j][i*mat->cols + k] = mat->array.SZ[j][k];
        }
     }
  }
  M->kgv = mat->kgv;
  M->flags.Integral = mat->flags.Integral;
  Check_mat(M);

  return(M);
}



/* -------------------------------------------------------------------- */
/* extract_r								*/
/* returns the spacegroup to the pointgroup G and the cozykel X		*/
/* code is taken out of extract.c                                       */
/* -------------------------------------------------------------------- */
bravais_TYP *extract_r(bravais_TYP *G,
                       matrix_TYP *X)
{
   bravais_TYP *H;

   int i, j;

   rat2kgv(X);
   Check_mat(X);

   /* is it a valid cocycle? */
   if ((G->dim * G->gen_no != X->rows) || (X->cols != 1)){
      fprintf(stderr,"The cocycle is not compatible to this point group\n");
      fprintf(stderr,"It should have %d * %d = %d rows\n",G->dim,G->gen_no,
                       G->dim*G->gen_no);
      exit(3);
   }

   H = init_bravais(G->dim+1);
   H->gen_no = G->gen_no;
   H->gen = (matrix_TYP **) malloc(G->gen_no * sizeof(matrix_TYP *));

   for (i=0;i<H->gen_no;i++){
      H->gen[i] = copy_mat(G->gen[i]);
      rat2kgv(H->gen[i]);
      Check_mat(H->gen[i]);
      real_mat(H->gen[i],H->dim,H->dim);
      iscal_mul(H->gen[i],X->kgv);
      H->gen[i]->kgv = H->gen[i]->kgv * X->kgv;
      for (j=0;j<H->dim-1;j++)
         H->gen[i]->array.SZ[j][H->dim-1] = X->array.SZ[i*(H->dim-1)+j][0];
      H->gen[i]->array.SZ[H->dim-1][H->dim-1] = X->kgv;
      Check_mat(H->gen[i]);
   }
  
   return(H);
}



/* -------------------------------------------------------------------- */
/* all_cocycles:                                                        */
/* Calculates the coycles for the spacegroup in the Z-class given by    */
/* the pointgroup G.                                                    */
/* -------------------------------------------------------------------- */
/* relator_input: presentation of G                                     */
/* G            : pointgroup representing the Z-class                   */
/*                CAUTION: the programm will only work correctly, if    */
/*                the correct normalizer of the pointgroup is given.    */
/* anzahl       : the function will set anzahl = number of affine       */
/*		  classes						*/
/* The following code is a part of the mainfunction of the programm     */
/* Extensions! 								*/
/* -------------------------------------------------------------------- */

matrix_TYP **all_cocycles(matrix_TYP *relator_input,
                          bravais_TYP *G,
                          int *anzahl,
                          matrix_TYP **matinv,
                          matrix_TYP ***X,
                          MP_INT **names,
                          int ****WORDS,
                          int **NUMBER_OF_WORDS,
                          matrix_TYP ***N,
                          int *coho_size,
                          int **list_of_names,
                          boolean l_option)
{
  matrix_TYP **Y;

  word *relator;

  int *len, i;

  long dim;

  MP_INT cohom_size;


/* we have to have at least the identity to generate the normalizer */
  if (G->normal_no == 0){
     G->normal_no = 1;
     G->normal = (matrix_TYP **) malloc(1 * sizeof(matrix_TYP *));
     G->normal[0] = init_mat(G->dim,G->dim,"1");
  }

/* speicher fuer die worte */
  relator = (word *) calloc(relator_input->rows,sizeof(word));

/* konvertieren der inputmatrix in relator-format */
  for (i=0;i<relator_input->rows;i++){
    matrix_2_word(relator_input,relator+i,i);
  }

  X[0] = cohomology(&dim,G->gen,matinv,relator,G->gen_no,relator_input->rows);

/* there is a special case to handle, which is the case that there
isn't a cohomology group at all */
  if (X[0][0]->cols <1){
     Y = (matrix_TYP **)malloc(sizeof(matrix_TYP *));
     Y[0] = init_mat(G->gen_no * G->dim,1,"");
     anzahl[0] = 1;
     N[0] = NULL;
     coho_size[0] = 1;
     names[0] = (MP_INT *) malloc(sizeof(MP_INT));
     mpz_init_set_si(names[0], 0);
  }
  else {
     cohom_size = cohomology_size(X[0][1]);
     coho_size[0] = mpz_get_ui(&cohom_size);

     if (l_option){
        list_of_names[0] = NULL;
     }
     else{
        if (coho_size[0] < TWOTO21){
           list_of_names[0] = (int *)calloc(coho_size[0], sizeof(int));
        }
        else{
           fprintf(stderr, "The cohomology group is too big!\n");
           fprintf(stderr, "If you are really interested in the graph of group - subgroup relations\n");
           fprintf(stderr, "please start the program again with the option -l!\n");
           exit(9);
        }
     }

     Y = extensions_o(X[0][0],X[0][1],X[0][2],G,&len,names,anzahl,
                      WORDS, NUMBER_OF_WORDS, N, cohom_size, 0, list_of_names[0]);
     free(len);
     mpz_clear(&cohom_size);
  }

  for (i=0;i<relator_input->rows;i++) wordfree(relator+i);
  free(relator);

  return(Y);
}



/* -------------------------------------------------------------------- */
/* transform a representative into the representation of                */
/* C_d1 x ... x C_dn                                                    */
/* -------------------------------------------------------------------- */
/* coz: cocycle in question                                             */
/* GLS: inverse of the third matrix returned by cohomology              */
/* D:   second matrix returned by cohomology                            */
/* -------------------------------------------------------------------- */
matrix_TYP *standard_rep(matrix_TYP *coz,
                         matrix_TYP *GLS,
                         matrix_TYP *D)
{
   int i, j,
       denominator,
       first, last, diff;

   matrix_TYP *tmp,
              *rep;


   for (first = 0; first < D->cols && D->array.SZ[first][first] == 1; first++);
   for (last = first; last < D->cols && D->array.SZ[last][last] != 0; last++);
   diff = last - first;
   tmp = mat_mul(GLS, coz);
   if (diff == 0)
      rep = init_mat(0, 0, "");
   else
      rep = init_mat(diff, 1, "");
   denominator = tmp->kgv;

   for (i = 0; i < diff; i++){
      j = tmp->array.SZ[i + first][0] * D->array.SZ[i + first][i + first];
      if ((j % denominator) != 0){
         fprintf(stderr,"ERROR in standard_rep: are you sure this is a cocycle?\n");
         fprintf(stderr,"If so, please report to the authors: carat@momo.math.rwth-aachen.de\n");
         exit(3);
      }
      else{
         rep->array.SZ[i][0] = (j / denominator) % D->array.SZ[i + first][i + first];
         while (rep->array.SZ[i][0] < 0)
            rep->array.SZ[i][0] += D->array.SZ[i + first][i + first];
      }
   }

   free_mat(tmp);

   return(rep);
}



/* -------------------------------------------------------------------- */
/* returns 1, if m in list; returns 0 otherwise                         */
/* -------------------------------------------------------------------- */
int yet_there(matrix_TYP *m,
              matrix_TYP **list,
              int no)
{
   int i;

   for (i = 0; i < no; i++)
      if (cmp_mat(m, list[i]) == 0)
         return(1);
   return(0);
}



/* -------------------------------------------------------------------- */
/* returns 1, if m = (0)_i,j, returns 0 otherwise                       */
/* -------------------------------------------------------------------- */
int equal_zero(matrix_TYP *m)
{
   int i, j;

   for (i = 0; i < m->rows; i++)
      for (j = 0; j < m->cols; j++)
         if (m->array.SZ[i][j] != 0)
            return(0);
   return(1);
}


/* -------------------------------------------------------------------- */
/* Generate list of matrices with the columns of M                      */
/* only the entries in M->array.SZ are copied                           */
/* -------------------------------------------------------------------- */
matrix_TYP **col_to_list(matrix_TYP *M)
{
   matrix_TYP **list;

   int i, j, rows, cols;


   cols = M->cols;
   rows = M->rows;
   list = (matrix_TYP **)calloc(cols, sizeof(matrix_TYP *));

   for (i = 0; i < cols; i++){
      list[i] = init_mat(rows, 1, "");
      for (j = 0; j < rows; j++){
         list[i]->array.SZ[j][0] = M->array.SZ[j][i];
      }
   }

   return(list);
}



/* -------------------------------------------------------------------- */
void free_H1_mod_ker_TYP(H1_mod_ker_TYP H1_mod_ker)
{
   int i;

   if (H1_mod_ker.flag == 0){
      for (i = 0; i < H1_mod_ker.erz_no; i++){
         free_mat(H1_mod_ker.M[i]);
      }
      free(H1_mod_ker.M);
      free_mat(H1_mod_ker.D);
      free_mat(H1_mod_ker.i);
   }
}



/* -------------------------------------------------------------------- */
static void mod(int *a, int b)
{
   a[0] = a[0] % b;
   if (a[0]<0) a[0] += b;
}



/* -------------------------------------------------------------------- */
/* return the number of the affine class for the cocylce coz            */
/* -------------------------------------------------------------------- */
/* data: Informationen ueber die Q-Klasse                               */
/* coz: Cozykel als Element von H^1                                     */
/* i: Nummer der Z-Klasse                                               */
/* flag: true => gebe -(Nummer der aff. Klasse) zurueck, falls nicht    */
/*       Standardvertreter der affinen Klasse                           */
/* wortflag: speichere ein Wort fuer Normalisatorelement, welches       */
/*           coz zu Standardverterter konjugiert                        */
/* wort: speicher Wort hier                                             */
/* -------------------------------------------------------------------- */
int number_of_affine_class(Q_data_TYP *data,
                           matrix_TYP *coz,
                           int i,
                           int flag,
			   boolean wortflag,
			   int **wort)
{
   matrix_TYP *rep;

   MP_INT nummer;

   int anz, n, first, last, zahl, pos;

   char *B;



   mpz_init(&nummer);
   for (first = 0; first < data->X[i][1]->cols && data->X[i][1]->array.SZ[first][first] == 1; first++);
   for (last = first; last < data->X[i][1]->cols && data->X[i][1]->array.SZ[last][last] != 0; last++);
   for (n = first; n < last; n++){
      mod(coz->array.SZ[n - first], data->X[i][1]->array.SZ[n][n]);
   }


   if (data->l_option == TRUE || wortflag){
      B = (char *)calloc(data->coho_size[i], sizeof(char));
      rep = orbit_rep(coz, data->N[i], data->Z[i]->normal_no, data->X[i][1], 0,
                      B, &nummer, &anz, wort, wortflag, NULL, NULL);

      pos = mpz_get_ui(&nummer);
      for (n = 0; n < data->aff_no[i]; n++){
         if (mpz_cmp(&nummer, &data->names[i][n]) == 0){
            break;
         }
      }
      if (n == data->aff_no[i]){
         fprintf(stderr, "ERROR in number_of_affine_class!\n");
         exit(6);
      }

      free_mat(rep);
      free(B);

      if (flag == 1){
         valuation(coz, data->X[i][1], &nummer);
         if (pos != mpz_get_ui(&nummer)){
            mpz_clear(&nummer);
            return(-n);
         }
      }
      mpz_clear(&nummer);
   }
   else{
      valuation(coz, data->X[i][1], &nummer);
      pos = mpz_get_ui(&nummer);
      zahl = data->list_of_names[i][pos];
      mpz_clear(&nummer);

      for (n = 0; n < data->aff_no[i]; n++){
         if (zahl == data->names_int[i][n]){
            break;
         }
      }
      if (n == data->aff_no[i]){
         fprintf(stderr, "ERROR in number_of_affine_class!\n");
         exit(6);
      }
      if (flag == 1 && pos != zahl)
         return(-n);
  }

  return(n);
}



/* -------------------------------------------------------------------- */
/* adds to all elements a row with "1"                                  */
/* -------------------------------------------------------------------- */
void kernel_elements_2_affine(matrix_TYP **elem,
                              int coho_size)
{
   int i, rows;


   rows = elem[0]->rows;

   for (i = 0; i < coho_size; i++){
      if (elem[i] != NULL){
         real_mat(elem[i], rows + 1, 1);
         elem[i]->array.SZ[rows][0] = 1;
      }
   }
}



/* -------------------------------------------------------------------- */
/* append G->cen to G->normal and set G->cen = NULL                     */
/* -------------------------------------------------------------------- */
void cen_to_norm(bravais_TYP *G)
{
   int i;

   if (G->cen_no > 0){
      G->normal = (matrix_TYP **)realloc(G->normal, 
                  (G->cen_no + G->normal_no) * sizeof(matrix_TYP *));
      for (i = 0; i < G->cen_no; i++){
         G->normal[G->normal_no + i] = G->cen[i];
      }
      free(G->cen);
      G->cen = NULL;
      G->normal_no += G->cen_no;
      G->cen_no = 0;
   }
}












