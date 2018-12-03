#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <tsubgroups.h>
#include <presentation.h>
#include <longtools.h>
#include <base.h>
#include <graph.h>
#include <bravais.h>
#include <datei.h>

int INFO_LEVEL;
extern int SFLAG;
boolean GRAPH_DEBUG;




int main (int argc, char *argv[])
{
   bravais_TYP *R;

   bravais_TYP **supergroups;

   int i, anz;

   char  comment[1000],
 	 file[1000];



   /* Daten und Optionen einlesen */
   /* =========================== */
   read_header(argc, argv);
   if (FILEANZ < 1 ||
       (is_option('h') && optionnumber('h') == 0)){
      printf("\n");
      printf("Usage: %s 'file1' [-h] [-f] [-t]\n", argv[0]);
      printf("\n");
      printf("Calculate the minimal translationengleich supergroups of a given\n");
      printf("space group up to conjugation under the affine normalizer of\n");
      printf("this space group.\n");
      printf("The program uses a database which only contains information\n");
      printf("up to dimension 4! So it will not work for higher dimensional space groups.\n");
      printf("\n");
      printf("file1: Spacegroup R in standard affine form (without translations)\n");
      printf("\n");
      printf("Options:\n");
      printf("-h : Give this help.\n");
      printf("-f : Write representatives in files `file1`_i.\n");
      printf("-t : Omit the last dim R generators.\n");
      printf("\n");
      printf("Cf.: TSubgroups\n\n");
      exit(11);
   }

   INFO_LEVEL = optionnumber('h');
   if (INFO_LEVEL & 12){
      SFLAG = 1;
   }

   /* get data */
   /* ======== */
   R = get_bravais(FILENAMES[0]);
   if (is_option('t')){
      R->gen_no -= (R->dim - 1);
   }

   supergroups = tsupergroups(R, &anz);
   for (i = 0; i < anz; i++){
      sprintf(comment, "%i-th supergroup of %s", i + 1, FILENAMES[0]);
      sprintf(file, "%s_%d", FILENAMES[0], i + 1);
      if (is_option('f'))
         put_bravais(supergroups[i], file, comment);
      else
         put_bravais(supergroups[i], 0, comment);
   }

   /* Speicherfreigabe */
   /* ================ */
   if (is_option('t')){
      R->gen_no += (R->dim - 1);
   }
   free_bravais(R);
   for (i = 0; i < anz; i++){
      free_bravais(supergroups[i]);
   }
   free(supergroups);

   /* Debugging */
   /* ========= */
   if (INFO_LEVEL & 12){
      fprintf(stderr,"write pointer_statistics\n");
      pointer_statistics(0,0);
   }

   exit(0);
}



