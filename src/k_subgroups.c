/* author: Oliver Heidbuechel */
/* last change: 07.02.2001 */

#include <ZZ.h>
#include<typedef.h>
#include<getput.h>
#include<matrix.h>
#include<longtools.h>
#include<tools.h>
#include"zass.h"
#include <base.h>
#include <bravais.h>
#include <graph.h>
#include <presentation.h>





int main (int argc, char *argv[])
{
   bravais_TYP *G, *R, **S;

   int i, anz, panz, OPT[6];

   char  comment[1000],
 	 file[1000];

   matrix_TYP **presentation,
              **base;

   bahn **strong;

 	
   read_header(argc, argv);
   if ( (is_option('p') && FILEANZ < 4) || (!is_option('p') && FILEANZ < 3) ||
        (is_option('h') && optionnumber('h') == 0) ){
      printf("\n");
      printf("Usage: %s 'file1' 'file2' [-p 'file3'] 'p_1' ... 'p_n' \n", argv[0]);
      printf("\n");
      printf("file1: Spacegroup R in standard affine form\n");
      printf("file2: Pointgroup G of R (G->gen and G->normal have to generate the normalizer)\n");
      printf("file3: (Optional) Presentation of G\n");
      printf("p_i  : primes\n");
      printf("\n");
      printf("Calculates the maximal klassengleich subgroups of R\n");
      printf("with p_i-power index.\n");
      printf("\n");
      printf("Options:\n");
      printf("-f   : print the subgroups in files 'file1_j'\n");
      printf("-h   : gives this help\n");
      printf("-d   : only for debugging, do not use\n");
      exit(11);
   }
   INFO_LEVEL = optionnumber('h');
   if (INFO_LEVEL & 12){
      SFLAG = 1;
   }

   /* get data */
   R = get_bravais(FILENAMES[0]);
   G = get_bravais(FILENAMES[1]);
   memset(G->divisors, 0, 100 * sizeof(int));
   if (is_option('p')){
      presentation = mget_mat(FILENAMES[2], &panz);
      if (panz > 1){
         fprintf(stderr, "you should only give a single matrix as presention\n");
         exit(3);
      }
      for (i = 3; i < FILEANZ; i++){
         G->divisors[atoi(FILENAMES[i])]++;
      }
   }
   else{
      base = get_base(G);
      strong = strong_generators(base,G,TRUE);
      presentation = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
      presentation[0] = pres(strong, G, OPT);
      for (i = 0; i < G->dim; i++){
         free_mat(base[i]);
         free_bahn(strong[i]);
         free(strong[i]);
      }
      free(strong);
      free(base);
      for (i = 2; i < FILEANZ; i++){
         G->divisors[atoi(FILENAMES[i])]++;
      }
   }

   /* calculate the minimal klassengleich subgroups of prime-power-index */
   S = max_k_sub(G, R, presentation[0], &anz, is_option('d'));

   /* print the minimal klassengleich subgroups of prime-power-index */
   if (is_option('f')){
      for (i = 0; i < anz; i++){
         sprintf(comment, "%d-th maximal klassengleich subgroup of %s", i + 1, FILENAMES[0]);
         sprintf(file, "%s_%d", FILENAMES[0], i + 1);
         put_bravais(S[i], file, comment);
      }
   }
   else{
      for (i = 0; i < anz; i++){
         sprintf(comment, "%d-th maximal klassengleich subgroup of %s", i + 1, FILENAMES[0]);
         put_bravais(S[i], 0, comment);
      }
   }

   /* clean */
   free_bravais(R);
   free_bravais(G);
   for (i = 0; i < anz; i++)
      free_bravais(S[i]);
   free(S);
   free_mat(presentation[0]);
   free(presentation);

   /* for debugging */
   if (INFO_LEVEL & 12){
      fprintf(stderr,"write pointer_statistics\n");
      pointer_statistics(0,0);
   }

   exit(0);
}






