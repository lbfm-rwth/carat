#include "typedef.h"
#include "gmp.h"
#include "getput.h"
#include "zass.h"
#include "matrix.h"
#include "tools.h"

void put_cocycle(matrix_TYP *COZ,
                 int dim,
                 int number,
                 char *file,
                 char *comment)
{
   int i,
       j;

   FILE *F;

   if (file == NULL){
      F = stdout;
   }
   else{
      F = fopen(file,"rw");
   }

   if (F == NULL){
      fprintf(stderr,"problems opening %s\n",file);
      exit(4);
   }

   if (dim * number != COZ->rows ){
      exit (3);
   }

   if (COZ->cols != 1){
      fprintf(stderr,"cozycle with more than 1 columns?\n");
      exit(3);
   }

   rat2kgv(COZ);
   Check_mat(COZ);

   if (COZ->kgv == 1 ||
       COZ->kgv == 0){
      fprintf(F,"%dx%d\t%s\n",dim,number,comment);
   }
   else{
      fprintf(F,"%dx%d/%d\t%s\n",dim,number,COZ->kgv,comment);
   }

   for (i=0;i<dim;i++){
      for (j=0;j<number;j++){
         fprintf(F,"%d ",COZ->array.SZ[j*dim+i][0]);
      }
      printf("\n");
   }

   return;
}

