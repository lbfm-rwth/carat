/* last change: 06.02.01 by Oliver Heidbuechel */


#include <ZZ.h>
#include <typedef.h>
#include <presentation.h>
#include <matrix.h>
#include <bravais.h>
#include <base.h>
#include <datei.h>
#include <graph.h>


/* -------------------------------------------------------------------- */
/* point_group								*/
/* extracts the linear part of the affine (space) group given in G.     */
/* Note this is the point group of the space group in case the space    */
/* group is given in standard form.   					*/
/* -------------------------------------------------------------------- */
bravais_TYP *p_group(bravais_TYP *G)
{
  bravais_TYP *H;

  int i;


  H = NULL;
  H = init_bravais(G->dim-1);

  H->gen_no = G->gen_no;
  H->gen = (matrix_TYP **) malloc(G->gen_no * sizeof(matrix_TYP *));
  for (i=0;i<H->gen_no;i++){
     H->gen[i] = copy_mat(G->gen[i]);
     real_mat(H->gen[i],H->dim,G->dim);
     real_mat(H->gen[i],H->dim,H->dim);
     rat2kgv(H->gen[i]);
     Check_mat(H->gen[i]);
  }
  return(H);
}



/* -------------------------------------------------------------------- */
/* plus_tranlationen:							*/
/* Adds translations to a spacegroup.					*/
/* -------------------------------------------------------------------- */
/* G  : the space group							*/
/* mat: matrix (the cols are the translating vectors)			*/
/* -------------------------------------------------------------------- */

void plus_translationen(bravais_TYP *G,
                        matrix_TYP *mat)
{
  matrix_TYP *hilfsmat;

  int i, j;


  /* create a matrix for help */
  hilfsmat = init_mat(G->dim, G->dim, "");
  for (i=0; i<hilfsmat->rows; i++){
     hilfsmat->array.SZ[i][i] = mat->kgv;
  }
  hilfsmat->kgv = mat->kgv;

  /* create the new translationmatrices for the spacegroup */
  G->gen = (matrix_TYP **)realloc(G->gen,(G->gen_no + mat->cols)*sizeof(matrix_TYP *));
  for (i=0; i<mat->cols; i++){
     G->gen[G->gen_no + i] = copy_mat(hilfsmat);
     for (j=0; j<G->dim-1; j++)
        G->gen[G->gen_no + i]->array.SZ[j][G->dim-1] =
           mat->array.SZ[j][i];
     Check_mat(G->gen[G->gen_no + i]);
  }
  G->gen_no = G->gen_no + mat->cols;

  /* clean up */
  free_mat(hilfsmat);
}



/* -------------------------------------------------------------------- */
/* extract_c								*/
/* extracts the translational part of the affine (space) group G as a   */
/* vector system (1-cocycle). 						*/
/* -------------------------------------------------------------------- */
matrix_TYP *extract_c(bravais_TYP *G)
{
  matrix_TYP *X;

  int i, kgv, j;

  X = NULL;


  X = init_mat(G->gen_no * (G->dim-1),1,"");
  kgv = 1;
  for (i=0;i<G->gen_no;i++){
     rat2kgv(G->gen[i]);
     kgv = kgv*G->gen[i]->kgv/GGT(kgv,G->gen[i]->kgv);
  }

  /* set the cocycle */
  X->kgv = kgv;
  for(i=0;i<G->gen_no;i++)
     for (j=0;j<G->dim-1;j++)
       X->array.SZ[i*(G->dim-1)+j][0] =  X->kgv/G->gen[i]->kgv *
                                       G->gen[i]->array.SZ[j][G->dim-1];
  Check_mat(X);
  return(X);
}











