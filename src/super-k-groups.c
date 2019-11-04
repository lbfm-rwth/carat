/* author: Oliver Heidbuechel */
/* last change: 09.12.2002 */

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
#include <name.h>





int main (int argc, char *argv[])
{
   bravais_TYP *P, *R, **S;

   int i, anz, OPT[6], *orbitlength = NULL;

   char  comment[1000],
 	 file[1000];

   bahn **strong;

   matrix_TYP *presentation,
              **base;


   read_header(argc, argv);
   if ( (is_option('p') && FILEANZ < 3) || (!is_option('p') && FILEANZ < 2) ||
        (is_option('h') && optionnumber('h') == 0) ){
      printf("\n");
      printf("Usage: %s 'file1' [-h] [-a [-n]] [-f] [-t] [-p 'file2'] 'p_1' ... 'p_n'\n", argv[0]);
      printf("\n");
      printf("file1: Spacegroup R in standard affine form (without translations)\n");
      printf("file2: (optional) Pointgroup P of R (with correct normalizer)\n");
      printf("p_i  : primes\n");
      printf("\n");
      printf("Calculates the minimal klassengleich supergroups with p_i-power index.\n");
      printf("\n");
      printf("Options:\n");
      printf("-h   : Give this help.\n");
      printf("-a   : Calculate supergroups up to conjugation under\n");
      printf("       the affine normalizer of R.\n");
      printf("-n   : Print the lengthes of the orbits\n");
      printf("       (this only works together with the option a).\n");
      printf("-f   : Print the supergroups in files 'file1_j'.\n");
      printf("-t   : Omit the last dim R generators.\n");
      printf("-d   : Only for debugging!\n");
      printf("\n");
      printf("Cf.: KSubgroups, Graph\n\n");
      exit(11);
   }
   INFO_LEVEL = optionnumber('h');
   if (INFO_LEVEL & 12){
      SFLAG = 1;
   }

   /* get data */
   R = get_bravais(FILENAMES[0]);
   if (is_option('t'))
      R->gen_no -= (R->dim - 1);
      
   if (is_option('p')){
      P = get_bravais(FILENAMES[1]);
   }
   else{
      P = point_group(R, 2);
      normalisator(P, NULL, NULL, 0, FALSE, FALSE);
   }

   /* presentation */
   base = get_base(P);
   strong = strong_generators(base, P, TRUE);
   OPT[0] = 0;
   presentation = pres(strong, P, OPT);
   for (i = 0; i < P->dim; i++){
      free_mat(base[i]);
      free_bahn(strong[i]);
      free(strong[i]);
   }
   free(strong);
   free(base);

   /* primes */
   memset(P->divisors, 0, 100 * sizeof(int));
   if (is_option('p')){
      for (i = 2; i < FILEANZ; i++)
         P->divisors[atoi(FILENAMES[i])]++;
   }
   else{
      for (i = 1; i < FILEANZ; i++)
         P->divisors[atoi(FILENAMES[i])]++;
   }

   /* calculate and print minimal klassengleich supergroups of prime-power-index */
   S = min_k_super(P, R, presentation, &anz, is_option('a'), is_option('d'), &orbitlength);

   /* print the minimal klassengleich subgroups of prime-power-index */
   if (is_option('a') && is_option('n')){
      fprintf(stderr, "\nLengthes of the orbits:\n");
      for (i = 1; i <= anz; i++)
         fprintf(stderr, "%i ", orbitlength[i]);
      fprintf(stderr, "\n\n");
   }
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
   if (orbitlength)
      free(orbitlength);
   if (is_option('t'))
      R->gen_no += (R->dim - 1);
   free_bravais(R);
   free_bravais(P);
   for (i = 0; i < anz; i++)
      free_bravais(S[i]);
   free(S);
   free_mat(presentation);

   /* for debugging */
   if (INFO_LEVEL & 12){
      fprintf(stderr,"write pointer_statistics\n");
      pointer_statistics(0,0);
   }

   exit(0);
}






