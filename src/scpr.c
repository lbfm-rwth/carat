#include "typedef.h"
#include "getput.h"
#include "matrix.h"


int main (int argc, char *argv[])

{  
   matrix_TYP **V,
              **F,
              **V_tr,
               *SCPR,
               *TT;

   char comment[1024];

   int truth,
       i,
       j,
       mv,
       nf;

   read_header(argc, argv);

   if (is_option('h') || FILEANZ ==0){
      printf("Usage: %s 'file1' ['file2']\n",argv[0]);
      printf("\n");
      printf("file1: matrix_TYP containing matrices (i. e. of basis transformations or group\n");
      printf("       elements)\n");
      printf("file2: (OPTIONAL) matrix_TYP containing Gram matrices\n");
      printf("\n");
      printf("Forms the scalar product of the COLUMNS of matrices A in file1 w.r.t. the\n");
      printf("Gram matrices F in file2, i. e. it computes A^(tr) * F * A.\n");
      printf("If file2 is ommitted, the identity matrix is used instead.\n");
      if (is_option('h')){
         exit(0);
      }
      else{
         exit(31);
      }
   }

   V = mget_mat (FILENAMES[0], &mv);

   V_tr = (matrix_TYP **) malloc(mv * sizeof(matrix_TYP *));

   for (j=0;j<mv;j++){
     V_tr[j] = tr_pose(V[j]);
   }

   if (FILEANZ > 1){
      F = mget_mat(FILENAMES[1],&nf);
   }
   else{
      nf = 1;
   }

   printf("#%d\n",mv*nf);
 
   for (i=0;i<nf;i++){
      for (j=0;j<mv;j++){
         if (FILEANZ>1){
            TT = mat_mul(F[i],V[j]);
            SCPR = mat_mul(V_tr[j],TT);
            free_mat(TT);
            sprintf(comment,"scalarproducts of %d-th matrix in %s w.r.t %s(cols)",
                              j+1,FILENAMES[0],FILENAMES[1]);
         }
         else{
            SCPR = mat_mul(V_tr[j],V[j]);
            sprintf(comment,"scalarproducts of %d-th matrix in %s(cols)",
                              j+1,FILENAMES[0]);
         }
         Check_mat(SCPR);
         put_mat(SCPR,NULL,comment,2);
         free_mat(SCPR);
      }
   }

   /* clean up memomry */
   if (FILEANZ > 1){
      for (i=0;i<nf;i++) free_mat(F[i]);
      free(F);
   }
   for (j=0;j<mv;j++){
      free_mat(V[j]);
      free_mat(V_tr[j]);
   }
   free(V);
   free(V_tr);

   exit(0);
}
