#include "typedef.h"
#include "bravais.h"
#include "matrix.h"


bravais_TYP *point_group(bravais_TYP *R,
                         int option)
{

   bravais_TYP *RES;

   int i;


   RES = init_bravais(R->dim-1);

   RES->gen_no = R->gen_no;

   RES->gen = (matrix_TYP **)  malloc(RES->gen_no * sizeof(matrix_TYP *));

   for (i=0;i<RES->gen_no;i++){
      RES->gen[i] = copy_mat(R->gen[i]);
      real_mat(RES->gen[i],RES->gen[i]->rows-1,RES->gen[i]->cols);
      real_mat(RES->gen[i],RES->gen[i]->rows,RES->gen[i]->cols-1);
      Check_mat(RES->gen[i]);
   }


   if (option & 2){
      /* calculate the space of invarinat forms, we might need it */
      RES->form = formspace(RES->gen,RES->gen_no,1,&RES->form_no);
   }

   return RES;

}


