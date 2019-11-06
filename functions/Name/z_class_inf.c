#include "ZZ.h"
#include "typedef.h"
#include "getput.h"
#include "name.h"
#include "gmp.h"
#include "name.h"
#include "bravais.h"
#include "datei.h"
#include "matrix.h"
#include "voronoi.h"
#include "autgrp.h"
#include "symm.h"
#include "contrib.h"
#include "base.h"
#include "zass.h"
#include "longtools.h"
#include "sort.h"
#include "idem.h"

static void rule_out_centerings(bravais_TYP *G,
                                matrix_TYP *L)
{

   int i;

   matrix_TYP *ele_L,
              *ele_tmp;

   /* ele_tmp = mat_inv(L); */
   ele_L = long_elt_mat(NULL,L,NULL);
   /* free_mat(ele_tmp); */

   for (i=0;i<G->zentr_no;i++){
      ele_tmp = long_elt_mat(NULL,G->zentr[i],NULL);
      if (mat_comp(ele_tmp, ele_L) != 0){
         free_mat(G->zentr[i]);
         G->zentr[i] = NULL;
      }
      free_mat(ele_tmp);
   }

   free_mat(ele_L);

   return;

}


matrix_TYP *z_class_inf(bravais_TYP *G,
                        bravais_TYP *DATABASEGROUP,
                        bravais_TYP **RES,
                        int *name)
{

   int i = 0,
       j,
       number,
       number2;

   bravais_TYP **QCLASS,
               **QCLASS2,
                *G_tr = NULL,
                *G_ALMOST;

   matrix_TYP *T = NULL,
              *ADLATTICE;

   /* split the Q-class of DATABASEGROUP into Z-classes, and
      then test each resulting one for Z-equivalence.
      Attention: the proper way to do this would be first to
      calculate only the homogenously decomposable groups,
      fits these, and then test invariant (i.e. index in the homogenously
      decomposable lattice), and save a lot of work. */
   DATABASEGROUP->form = formspace(DATABASEGROUP->gen,
                                   DATABASEGROUP->gen_no,
                                   1,
                                   &DATABASEGROUP->form_no);
   QCLASS = q2z(DATABASEGROUP,&number,TRUE, NULL, TRUE, FALSE);

   /* calculate the almost decomposable lattice */
   ADLATTICE = almost_decomposable_lattice(G);
   G_ALMOST = konj_bravais(G,ADLATTICE);

   /* added tilman 28.02.00 */
   if (G_ALMOST->form == NULL){
      G_ALMOST->form = formspace(G_ALMOST->gen,
                                 G_ALMOST->gen_no,1,&G_ALMOST->form_no);
   }

   long_rein_formspace(G_ALMOST->form,G_ALMOST->form_no,1);

   if (INFO_LEVEL & 4){
      put_mat(ADLATTICE,0,0,0);
   }

   /* G_ALMOST is a proper group, ie. has the space of invariant
      forms assigned */
   while(T == NULL && i < number){
      T = z_equivalent(G_ALMOST, &G_tr , QCLASS[i]);

      if (T){
         free_bravais(G_tr); G_tr = NULL;
         free_mat(T); T = NULL;

         rule_out_centerings(QCLASS[i],ADLATTICE);
         QCLASS2 = get_groups(QCLASS+i,1,&number2);
         j = 0;
         while(T == NULL && j < number2){
            if (QCLASS2[j]){
               T = z_equivalent(G, &G_tr , QCLASS2[j]);
            }
            j++;
         }
         if (T == NULL){
            fprintf(stderr,"error in z_class_inf.\n");
            fprintf(stderr,"group does not appear in catalog.\n");
            put_bravais(G,0,0);
            exit(4);
         }
      }
      i++;
   }
   free_bravais(G_tr);
   free_bravais(G_ALMOST);

   if (T == NULL){
      fprintf(stderr,"error in z_class_inf.\n");
      fprintf(stderr,"group does not appear in catalog.\n");
      put_bravais(G,0,0);
      exit(4);
   }
   else{
      name[0] = i;
      name[1] = j;
      *RES = QCLASS2[j-1];
   }

   for (i=0;i<number;i++){
      free_bravais(QCLASS[i]);
   }
   for (i=0;i<number2;i++){
      if (QCLASS2[i] != NULL && *RES != QCLASS2[i]) free_bravais(QCLASS2[i]);
   }
   free(QCLASS);
   free(QCLASS2);
   free_mat(ADLATTICE);

   return T;
}


