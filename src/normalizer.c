#include "typedef.h"
#include "bravais.h"
#include "base.h"
#include "matrix.h"
#include "getput.h"
#include "datei.h"

int INFO_LEVEL;
extern int SFLAG;

int main (int argc, char *argv[])
{
   bravais_TYP *Gtr,
               *H;

   matrix_TYP *A;
   
   int prime;


   read_header(argc, argv);

   if((FILEANZ < 1) || (is_option('h') && optionnumber('h')==0))
   {
     printf("Usage: Normalizer 'file1' ['file2' ['file3']] [-p=prime] [-o]\n");
     printf("\n");
     printf("file1: bravais_TYP containing the group G.\n");
     printf("file2: (OPTIONAL) bravais_TYP containing the transposed of G. (cf. Tr_bravais)\n");
     printf("file3: (OPTIONAL) matrix_TYP of a G-perfect form. (cf. First_perfect)\n");
     printf("\n");
     printf("Calculates a set of matrices which together with G generate the\n");
     printf("normalizer N_GL_n(Z) (G) of G in GL_n(Z).\n");
     printf("NOTE: the output echoes any input information about G, except input about\n");
     printf("generators of the normalizer.\n");
     printf("NOTE: The dimension of the space of invariant forms is a measure for the\n");
     printf("complexity of the algorithm. Up to degree 6 the only infeasible case are\n");
     printf("<I_6> and <-I_6>. Here the generators of the normalizer can be taken\n");
     printf("from `Bravais_cat' with family 1,1,1,1,1,1.\n");
     printf("\n");
     printf("Options:\n");
     printf("-b      :  The normalizer of the bravais group B(G) is calculated. With this\n");
     printf("           option the program is much faster. (The normalizer of G is a\n");
     printf("           subgroup of N_GL_n(Z) (B(G)). )\n");
     printf("-p=prime:  The determinants of the perfect forms are\n");
     printf("           calculated module prime. The default is 1949.\n");
     printf("-o      :  The G-perfect forms are given as additional output.\n");
     printf("\n");
     if (is_option('h')){
        exit(0);
     }
     else{
        exit(31);
     }
   }

   INFO_LEVEL = optionnumber('h');
   if (INFO_LEVEL & 12){
      SFLAG = 1;
   }


   /* get data */
   H = get_bravais(FILENAMES[0]);

   if (FILEANZ > 1){
      Gtr = get_bravais(FILENAMES[1]);
   }
   else{
      Gtr = NULL;
   }

   /* read an G-perfect form if it is given */
   if (FILEANZ > 2){
      A = get_mat(FILENAMES[2]);
   }
   else{
      A = NULL;
   }

   prime = optionnumber('p');


   /* calculate normalizer */
   normalisator(H, Gtr, A, prime, is_option('b'), is_option('o'));
   if(!is_option('o')){
      put_bravais(H, NULL, "group with complete normalizer");
   }


   /* cleaning up the memory */
   if (A != NULL)
      free_mat(A);
   free_bravais(H);
   if (Gtr != NULL)
      free_bravais(Gtr);

   /* some diagnostic for memory leakage */
   if (INFO_LEVEL & 12){
      pointer_statistics(0,0);
   }

   exit(0);
}

