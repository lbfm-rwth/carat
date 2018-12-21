#include "ZZ.h"
#include "typedef.h"
#include "getput.h"
#include "gmp.h"
#include "name.h"
#include "bravais.h"
#include "datei.h"
#include "matrix.h"
#include "voronoi.h"
#include "autgrp.h"
#include "symm.h"
#include "base.h"
#include "zass.h"
#include "longtools.h"



void extend(matrix_TYP *T)
{

   real_mat(T,T->rows,T->cols+1);
   real_mat(T,T->rows+1,T->cols);

   T->array.SZ[T->rows-1][T->cols-1] = T->kgv;

   return;
}


static matrix_TYP *get_cocycle(bravais_TYP *R,
                               bravais_TYP *P)
{


   matrix_TYP **RG,
               *coz;

   int  k,
        j,
        denominator,
      **words;

   RG = (matrix_TYP **) malloc( R->gen_no * sizeof(matrix_TYP*));
   words = (int **) malloc(P->gen_no * sizeof(int *));
   denominator = 1;
   for (j=0;j<R->gen_no;j++){
       RG[j] = copy_mat(R->gen[j]);
       RG[j]->cols--;
       RG[j]->rows--;
       Check_mat(RG[j]);
       if (!RG[j]->flags.Integral){
          fprintf(stderr,"The point group has to be integral\n");
          exit(3);
       }
       rat2kgv(R->gen[j]);
       denominator *= (R->gen[j]->kgv / GGT(R->gen[j]->kgv,denominator));
     }

     /* stick the rigth INTEGRAL cozycle at the end of the RG[j] */
     for (j=0;j<R->gen_no;j++){
        RG[j]->cols++;
        RG[j]->rows++;
        for (k=0;k<RG[j]->rows-1;k++)
           RG[j]->array.SZ[k][R->dim-1] = (denominator / R->gen[j]->kgv) *
                   R->gen[j]->array.SZ[k][R->dim-1];
        RG[j]->array.SZ[R->dim-1][R->dim-1] = 1;
        Check_mat(RG[j]);
     }

     /* get the cozycle on the right generators */
     coz = reget_gen(RG,R->gen_no,P,words,TRUE);

     /* the cozykle has to become the right denominator */
     coz->kgv = denominator;
     Check_mat(coz);

     for (j=0;j<R->gen_no;j++){
        free_mat(RG[j]);
     }
     for (j=0;j<P->gen_no;j++){
        free(words[j]);
     }
     free(words);
     free(RG);

     return coz;

}


bravais_TYP *space_group_from_matrix(bravais_TYP *G,
                                     matrix_TYP *x,
                                     matrix_TYP *cocycle,
                                     matrix_TYP *D)
{

   bravais_TYP *R;

   matrix_TYP *C;

   int i,
       j,
       k;

   R = init_bravais(G->dim+1);
   C = convert_to_cozycle(x,cocycle,D);

   R->gen = (matrix_TYP **) malloc(G->gen_no * sizeof(matrix_TYP *));
   R->gen_no = G->gen_no;
   for (i=0;i<G->gen_no;i++){
      R->gen[i] = copy_mat(G->gen[i]);
      real_mat(R->gen[i],G->dim+1,G->dim+1);
      R->gen[i]->array.SZ[G->dim][G->dim] = 1;
      iscal_mul(R->gen[i],C->kgv);
      R->gen[i]->kgv = C->kgv;
   }

   /* stick the cocycle in the last column of the matrices generating R */
   k = 0;
   for (i=0;i<R->gen_no;i++){
      for (j=0;j<G->dim;j++){
         R->gen[i]->array.SZ[j][G->dim] = C->array.SZ[k][0];
         k++;
      }
      Check_mat(R->gen[i]);
   }

   free_mat(C);

   return R;

}

matrix_TYP *aff_class_inf(bravais_TYP *R,
                          bravais_TYP *DATAZ,
                          matrix_TYP *PRES,
                          MP_INT *aff_name,
                          bravais_TYP **RC)
{

   matrix_TYP *cozycle,
             **X,
             **Y,
              *coz_mat,
             **matinv,
              *RES;

   word *relator;

   int i;

   long dim;

   /* first thing to do is to find the generators of DATAZ in the
      point group of R */
   cozycle = get_cocycle(R,DATAZ);

   /* do the cohomology calculations */
   relator = (word *) calloc(PRES->rows,sizeof(word));
   for (i=0;i<PRES->rows;i++){
     matrix_2_word(PRES,relator+i,i);
   }


   matinv = (matrix_TYP **) calloc(DATAZ->gen_no , sizeof(matrix_TYP *));
   X = cohomology(&dim,DATAZ->gen,matinv,relator,DATAZ->gen_no,PRES->rows);

   if (X[0]->cols > 0){ 
      /* give the group a name */
      Y = identify(X[0],X[1],X[2],DATAZ,&cozycle,aff_name,1,3,NULL,NULL);
      RES=Y[0]; free(Y);

      if (RC){
         /* construct our representative */
         coz_mat = reverse_valuation(aff_name,X[1]);
         *RC = space_group_from_matrix(DATAZ,coz_mat,X[0],X[1]);
         free_mat(coz_mat);
      }

   }
   else{
      mpz_set_si(aff_name,0);
      RES = init_mat(DATAZ->dim+1,DATAZ->dim+1,"1");
      coboundary(DATAZ,cozycle,RES);

      if (RC){
         /* construct the split extension */
         RC[0] = init_bravais(DATAZ->dim + 1);
         RC[0]->gen = (matrix_TYP **) malloc(DATAZ->gen_no*sizeof(matrix_TYP*));
         RC[0]->gen_no = DATAZ->gen_no;
         for (i=0;i<DATAZ->gen_no;i++){
            RC[0]->gen[i] = copy_mat(DATAZ->gen[i]);
            extend(RC[0]->gen[i]);
         }
         free_mat(coz_mat);
      }
   }

   /* clean up and return */
   for (i=0;i<3;i++) free_mat(X[i]); free(X);
   free_mat(cozycle);
   for (i=0;i<PRES->rows;i++) wordfree(relator+i);
   free(relator);
   for (i=0;i<DATAZ->gen_no;i++)
      if (matinv[i] != NULL) free_mat(matinv[i]);
   free(matinv);

   return RES;
}


