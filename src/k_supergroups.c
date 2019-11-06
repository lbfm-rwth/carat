/* author: Oliver Heidbuechel */
/* last change: 19.01.2001 */

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
   bravais_TYP *P, *R, **S;

   int i, anz;

   char  comment[1000],
 	 file[1000];

 	
   read_header(argc, argv);
   if (FILEANZ < 1 || (is_option('h') && optionnumber('h') == 0)){
      printf("\n");
      printf("Usage: %s 'file1' 'p_1' ... 'p_n' \n", argv[0]);
      printf("\n");
      printf("file1: Spacegroup in standard affine form\n");
      printf("p_i  : primes\n");
      printf("\n");
      printf("Calculates the minimal klassengleich supergroups with p_i-power index.\n");
      printf("\n");
      printf("Options:\n");
      printf("-t   : omit the last dim generators of the spacegroup\n");
      printf("-f   : print the supergroups in files 'file1_j'\n");
      printf("-h   : gives this help\n");
      printf("-d   : only for debugging\n");
      exit(11);
   }
   INFO_LEVEL = optionnumber('h');

   /* get data */
   R = get_bravais(FILENAMES[0]);
   if (is_option('t'))
      R->gen_no -= (R->dim - 1);
   P = p_group(R);
   for (i = 2; i <= FILEANZ; i++){
      P->divisors[atoi(argv[i])]++;
   }

   /* calculate and print minimal klassengleich supergroups of prime-power-index */
   S = min_k_super(P, R, &anz, is_option('d'));

   if (is_option('f')){
      for (i = 0; i < anz; i++){
         sprintf(comment, "%d-th minimal klassengleich supergroup of %s", i + 1, FILENAMES[0]);
         sprintf(file, "%s_%d", FILENAMES[0], i + 1);
         put_bravais(S[i], file, comment);
      }
   }
   else
      for (i = 0; i < anz; i++)
         put_bravais(S[i],0,0);

   /* clean */
   if (is_option('t'))
      R->gen_no += (R->dim - 1);
   free_bravais(R);
   free_bravais(P);
   for (i = 0; i < anz; i++)
      free_bravais(S[i]);
   free(S);

   exit(0);
}






