#include "ZZ.h"
#include "typedef.h"
#include "voronoi.h"
#include "symm.h"
#include "autgrp.h"
#include "bravais.h"
#include "base.h"
#include "idem.h"
#include "longtools.h"
#include "reduction.h"


bravais_TYP **get_groups(bravais_TYP **ADGROUPS,
                         int ad_no,
                         int *number)
{
   int i,
       j,
       k,
       normal_no,
       centr_no;

   bravais_TYP **GROUPS,
                *G1,
                *G2;

   matrix_TYP  *tmp,
              **centerings,
              **normal;

   G1 = init_bravais(ADGROUPS[0]->dim);

   /* look for all the centerings and so on */
   number[0] = 0;
   for (i=0;i<ad_no;i++){
     number[0] += ADGROUPS[i]->zentr_no;
   }
   GROUPS = (bravais_TYP **) calloc(2 * number[0] , sizeof(bravais_TYP *));

   k=0;
   for (i=0;i<ad_no;i++){
      centerings = ADGROUPS[i]->zentr;
      centr_no = ADGROUPS[i]->zentr_no;
      ADGROUPS[i]->zentr = NULL;
      ADGROUPS[i]->zentr_no = 0;
      normal = ADGROUPS[i]->normal;
      ADGROUPS[i]->normal = NULL;
      normal_no = ADGROUPS[i]->normal_no;
      ADGROUPS[i]->normal_no = 0;
      G1->gen = normal;
      G1->gen_no = normal_no;
      for (j=0;j<centr_no;j++){
         if (centerings[j]){
            tmp = mat_inv(centerings[j]);
            if (j==0){
               GROUPS[k] = konj_bravais(ADGROUPS[i],tmp);
               GROUPS[k]->normal = normal;
               GROUPS[k]->normal_no = normal_no;
               GROUPS[k + number[0] ] = (bravais_TYP * ) 1;
            }
            else{
               G2 = gittstab(G1,centerings[j]);
               ADGROUPS[i]->normal = G2->gen;
               ADGROUPS[i]->normal_no = G2->gen_no;
               GROUPS[k] = konj_bravais(ADGROUPS[i],tmp);
               long_rein_formspace(GROUPS[k]->form,GROUPS[k]->form_no,1);
               ADGROUPS[i]->normal = NULL;
               ADGROUPS[i]->normal_no = 0;
               free_bravais(G2);
               GROUPS[k + number[0] ] = (bravais_TYP * ) 0;
            }
            free_mat(tmp);
            free_mat(centerings[j]);
         }
         k++;
      }
      if (centerings[0] == NULL){
         /* case only used for z_class_inf */
         ADGROUPS[0]->normal = normal;
         ADGROUPS[0]->normal_no = normal_no;
      }
      free(centerings);
   }

   free(G1);

   return GROUPS;
}
