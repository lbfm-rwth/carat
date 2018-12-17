#include "typedef.h"
#include "longtools.h"
#include "sort.h"
#include "matrix.h"
#include "orbit.h"
#include "idem.h"
#include "datei.h"
#include "tools.h"
#include "getput.h"

extern int INFO_LEVEL;

/******************************************************************************
@
@------------------------------------------------------------------------------
@ FILE: v4_catalog.c
@------------------------------------------------------------------------------
@
*******************************************************************************/

/******************************************************************************
@
@------------------------------------------------------------------------------
@
@ static matrix_TYP *test_v4(bravais_TYP *G)
@
@ Tests whether the given Bravais group G is isomorphic to V_4 = C_2 x C_2.
@ The order of the group has to be 4, and this order has to sit in G->order.
@
@ The return is NULL if the group is not isomophic to V_4, and otherwise
@ a matrix with non negative trace in G.
@------------------------------------------------------------------------------
@
*******************************************************************************/
static matrix_TYP *test_v4(bravais_TYP *G)
{
   static int i,
           j,
           opt[6];

   matrix_TYP **ELE,
               *ID,
               *TMP,
               *g;

   if (G->order != 4){
      return NULL;
   }

   opt[2] = 4;
   j = 0;

   /* calculate all elements */
   ID = init_mat(G->dim,G->dim,"1");
   ELE = orbit_alg(ID,G,NULL,opt,&i);

   for (i=0;i<4;i++){
      TMP = mat_mul(ELE[i],ELE[i]);
      if (mat_comp(ID,TMP)==0)
         j++;
      free_mat(TMP);
   }

   g = NULL;

   /* all elements should have order at most 2 */
   if (j==4){
      for (i=0;i<4 && g == NULL; i++){
         j = trace(ELE[i]);
         if ((0<=j) && (j <G->dim)){
            g = ELE[i];
            ELE[i] = NULL;
         }
      }
   }

   /* clean up */
   for (i=0;i<4;i++){
     if (ELE[i] != NULL) free_mat(ELE[i]);
   }
   free(ELE);

   free_mat(ID);

   return g;
}

/******************************************************************************
@
@------------------------------------------------------------------------------
@
@ static matrix_TYP *basis_part(matrix_TYP *g, matrix_TYP *h,
@                               matrix_TYP **plse,int rank2,int eps)
@
@------------------------------------------------------------------------------
@
*******************************************************************************/
static matrix_TYP *basis_part(matrix_TYP *g, matrix_TYP *h,
                              matrix_TYP **plse,int rank2,int eps)
{

  int i,
      j,
      k,
      old_rank,
      rank;

  matrix_TYP *E,
             *A,
            **tmp,
             *tmp2,
             *zass_mat,
             *RES;

  /* calculate a Z-basis for the eigen space of 0 */
  tmp = long_solve_mat(h,NULL);
  E = tmp[1];
  free(tmp);

  /* avoid trouble if there isn't any psle.. */
  if (rank2 == 0){
      return E;
  }

  /* calculate a Z-basis for the intersection E ^ (e1 + (eps) * g e1 .....) */
  zass_mat = init_mat(rank2,2*g->cols,"");
  for (i=0;i<rank2;i++){
     tmp2 = mat_mul(g,plse[i]);
     for (j=0;j<g->cols;j++){
         zass_mat->array.SZ[i][j] = plse[i]->array.SZ[j][0]
                                    + eps * tmp2->array.SZ[j][0];
         zass_mat->array.SZ[i][j+g->cols] = zass_mat->array.SZ[i][j];
     }
     free_mat(tmp2);
  }
  tmp2 = long_rein_mat(zass_mat);
  free_mat(zass_mat);
  zass_mat = tmp2;
  tmp2 = NULL;
  old_rank = zass_mat->rows;
  real_mat(zass_mat,zass_mat->rows+E->cols,zass_mat->cols);
  for (i=0;i<E->cols;i++){
     for (j=0;j<E->rows;j++){
         zass_mat->array.SZ[old_rank+i][j] = E->array.SZ[j][i];
     }
  }
  long_row_hnf(zass_mat);

  rank = -1;
  i = -1;
  while (rank == -1){
     i++;
     rank = i;
     for (j=0;j<g->cols;j++)
        if (zass_mat->array.SZ[i][j] != 0){
           rank = -1;
        }
  }

  tmp2 = init_mat(zass_mat->rows-rank,g->cols,"");
  for (i=0;i<zass_mat->rows-rank;i++){
     for (j=0;j<g->cols;j++){
        tmp2->array.SZ[i][j] = zass_mat->array.SZ[i+rank][j+g->cols];
     }
  }
  free_mat(zass_mat);
  zass_mat = long_rein_mat(tmp2);
  free_mat(tmp2);
  tmp2 = zass_mat;
  zass_mat=NULL;

  /* now add vectors to give a Z-basis of E */
  /* represent the rows of tmp2 as linear combinations of the columns of E */
  /* write these in the rows of A */
  zass_mat = tr_pose(tmp2);
  free_mat(tmp2);
  tmp = long_solve_mat(E,zass_mat);
  if (tmp[1] != NULL){
      fprintf(stderr,"error in basis_part\n");
      exit(3);
  }
  A = tr_pose(tmp[0]);
  free_mat(tmp[0]);
  free(tmp);
  free_mat(zass_mat);

  /* make an gauss as far as possible and add indices */
  old_rank = A->rows;
  RES = init_mat(g->cols,E->cols-old_rank,"");
  for (i=old_rank;i<E->cols;i++){
     long_row_hnf(A);
     real_mat(A,A->rows+1,A->cols);

     /* search for a "step" */
     for (j=0;j<A->rows && A->array.SZ[j][j]!=0;j++);
     j--;

     /* add one row vector of A and make sure the lattice stays pure */
     if (j < 0){
        A->array.SZ[A->rows-1][0] = 1;
     }
     else{
        gcd_darstell(A->array.SZ[j][j],A->array.SZ[j][j+1],
                  A->array.SZ[A->rows-1]+j+1,A->array.SZ[A->rows-1]+j,&k);
        A->array.SZ[A->rows-1][j+1] *= (-1);
     }

     /* stick this linear combination into RES */
     for (j=0;j<A->cols;j++){
        for (k=0;k<RES->rows;k++){
           RES->array.SZ[k][i-old_rank] += A->array.SZ[A->rows-1][j]
                                         * E->array.SZ[k][j];
        }
     }
  }

  /* cleanup */
  free_mat(A);
  free_mat(E);

  return RES;
}

/******************************************************************************
@
@------------------------------------------------------------------------------
@
@ static matrix_TYP *normalize(matrix_TYP *g,int *dim,int flag)
@
@------------------------------------------------------------------------------
@
*******************************************************************************/
static matrix_TYP *normalize(matrix_TYP *g,int *dim,int flag)
{

  matrix_TYP *h,
             *RES,
            **plse,
             *tmp;

  int i,
      j,
      k,
      rank2;

  /* set h = g - 1 */
  h = copy_mat(g);
  for (i=0;i<h->cols;i++)
    h->array.SZ[i][i]--;

  tmp = copy_mat(h);
  modp_mat(tmp,2);
  init_prime(2);
  rank2 = p_gauss(tmp);
  cleanup_prime();

  plse = (matrix_TYP **) malloc(rank2 * sizeof(matrix_TYP *));
  for (i=0;i<rank2;i++){
     plse[i] = init_mat(h->rows,1,"");
     k = 0;
     for (j=0;j<tmp->cols && k==0;j++)
        k = plse[i]->array.SZ[j][0] = tmp->array.SZ[i][j];
  }
  free_mat(tmp);

  if (flag || rank2 == dim[0]){
     dim[0] = rank2;

     RES = init_mat(g->cols,g->cols,"");

     for (i=0;i<rank2;i++){
        for (j=0;j<g->cols;j++){
           RES->array.SZ[j][2*i] = plse[i]->array.SZ[j][0];
        }
        tmp = mat_mul(g,plse[i]);
        for (j=0;j<g->cols;j++){
           RES->array.SZ[j][2*i+1] = tmp->array.SZ[j][0];
        }
        free_mat(tmp);
     }

     k = 2*rank2;
     tmp = basis_part(g,h,plse,rank2,1);
     for (i=0;i<tmp->cols;i++){
        for (j=0;j<tmp->rows;j++){
           RES->array.SZ[j][k] = tmp->array.SZ[j][i];
        }
        k++;
     }
     free_mat(tmp);

     for (i=0;i<h->cols;i++)
        h->array.SZ[i][i] += 2;

     tmp = basis_part(g,h,plse,rank2,-1);
     for (i=0;i<tmp->cols;i++){
        for (j=0;j<tmp->rows;j++){
           RES->array.SZ[j][k] = tmp->array.SZ[j][i];
        }
        k++;
     }
     free_mat(tmp);
  }
  else{
     /* we didn't find an appropriate element */
     RES = NULL;
  }

  /* free the memory */
  for (i=0;i<rank2;i++)
     free_mat(plse[i]);
  if (plse != NULL) free(plse);
  free_mat(h);

  if (RES != NULL){
     Check_mat(RES);
  }
 
  return RES;
}

/******************************************************************************
@
@------------------------------------------------------------------------------
@
@ bravais_TYP *catalog_number_v4(bravais_TYP *G,char *symb,
@                             matrix_TYP **TR,int *almost,int *zclass)
@
@ Tests wheter the Bravais group G is isomorphic to C_2 or V_4.
@ In these cases it will return a group which is Z-equivalent to G
@ from the catalog which sits in TOPDIR/tables.
@
@ It will return a transfromation matrix via TR[0], and the position
@ in the catalog via almost[0] and zclass[0].
@
@ bravais_TYP *G : The group in question. Its order must be given.
@ char *symb     : The symbol of the group. It can be calculated via symbol(..)
@ matrix_TYP **TR: pointer for the transformation matrix which transforms
@                  the given group G to the group returned via konj_bravais,
@                  ie. TR[0]  * G * TR[0]^-1 = group returned.
@ int *almost    : the position of the almost decomposable group in the
@                  catalog is returned via this pointer.
@ int *zclass    : 2 coordinate of the group in the catalog.
@
@------------------------------------------------------------------------------
@
*******************************************************************************/
bravais_TYP *catalog_number_v4(bravais_TYP *G,char *symb,
                            matrix_TYP **TR,int *almost,int *zclass)
{

   bravais_TYP *T = NULL,
               *H;

   symbol_out *S;

   matrix_TYP *X,
              *g,
              *g2,
              *h,
              *h2,
              *TG,
              *TG2;

   char *file;

   int i,
       dim_g_F2;

   if (G->dim > MAXDIM){
      fprintf(stderr,"This program does only work up to dimension %d\n",
                      MAXDIM);
      exit(3);
   }

   /* check whether we do have the trivial case of <-I_n> */
   if ((strcmp(symb,"1") == 0) ||
       (strcmp(symb,"1,1") == 0) || 
       (strcmp(symb,"1,1,1") == 0) || 
       (strcmp(symb,"1,1,1,1") == 0) || 
       (strcmp(symb,"1,1,1,1,1") == 0) || 
       (strcmp(symb,"1,1,1,1,1,1") == 0)){
       almost[0] = 1;
       zclass[0] = 1;
       S = read_symbol_from_string(symb);
       T = S->grp;
       free(S);
       free(S->fn);
       TR[0] = init_mat(G->dim,G->dim,"1");
       if (T->zentr_no > 0){
          fprintf(stderr,"catalog_number: an error ocurred in the -I_n case\n");
          exit(3);
       }
       return T;
   }

   /* firstly test whether we have a group isomorphic to v4 */
   g = test_v4(G);

   if (g!= NULL){
      /* we got an element of G of order 2, which has non negative trace */

      /* we should only have 1 almost decomposable family */
      S = read_symbol_from_string(symb);
      H = S->grp;
      get_zentr(S);
      if (S->fn != NULL){
         fprintf(stderr,"catalog_number: an error occured in v4-part\n");
         exit(3);
      }
      free(S);

      /* throw away centralizer & normalizer, ist only hinders */
      if (H->cen_no >0){
         for (i=0;i<H->cen_no;i++) free_mat(H->cen[i]);
         free(H->cen);
         H->cen = NULL;
         H->cen_no = 0;
      }
      if (H->normal_no >0){
         for (i=0;i<H->normal_no;i++) free_mat(H->normal[i]);
         free(H->normal);
         H->normal = NULL;
         H->normal_no = 0;
      }

      /* normalize g */
      TG = normalize(g,&dim_g_F2,TRUE);
      X = long_mat_inv(TG);
      TG2 = mat_kon(X,g,TG);
      free_mat(g);
      free_mat(X);
      g = TG2;

      if (INFO_LEVEL & 4){
         put_mat(TG,NULL,"TG",2);
         put_mat(g,NULL,"g",2);
      }

      /* bare in mind that the trace is Q-invariant */
      g2 = test_v4(H);

      almost[0] = 1;
      zclass[0] = -1;
      TG2 = NULL;
      while (TG2 ==0){

         if (zclass[0] == (-1)){
            h = copy_mat(g2);
         }
         else if (zclass[0] >= H->zentr_no){
            fprintf(stderr,"bravais_catalog: v4-part: didn't find the group\n");
            exit(3);
         }
         else{
            X = mat_inv(H->zentr[zclass[0]]);
            h = mat_kon(X,g2,H->zentr[zclass[0]]);
            free_mat(X);
         }

         TG2 = normalize(h,&dim_g_F2,FALSE);

         if (TG2 != NULL){
            /* we found the corresponding group */
            /* the transforming matrix is the quotient of TG and TG2 */
            if (INFO_LEVEL & 4) put_mat(TG2,NULL,"TG2",2);
            X = long_mat_inv(TG);
            TR[0] = mat_mul(TG2,X);
            if (zclass[0] < 0){
               T = H;
               H = NULL;
            }
            else{
               T = Z_class(H,H->zentr[zclass[0]]);
            }
         }

         free_mat(h);

         zclass[0]++;
      }

      /* we started off with -1 */
      zclass[0]++;

      free_mat(g);
      free_mat(g2);
      free_mat(TG);
      free_mat(TG2);
      free_mat(X);
      if (H != NULL) free_bravais(H);

   }

   return T;
}
