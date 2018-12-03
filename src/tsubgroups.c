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
   bravais_TYP *R,
               *P;

   TSubgroup_TYP **subgroups;

   matrix_TYP **presentation,
              **base,
              **gapwords;

   char  comment[1000],
 	 file[1000];

   bahn **strong;

   int OPT[6], panz, no, i, no_neu;


   /* Daten und Optionen einlesen */
   /* =========================== */
   read_header(argc, argv);
   if ((!is_option('g') && FILEANZ < 1) ||
       (is_option('g') && FILEANZ < 3) ||
       (is_option('h') && optionnumber('h') == 0) ){
      printf("\n");
      printf("Usage:\n");
      printf("1) %s 'file1' [-h] [-a] [-f] [-t]\n", argv[0]);
      printf("2) %s -g [-c] 'file1' 'file2' 'file3' ['file4']\n", argv[0]);
      printf("\n");
      printf("Calculate the maximal translationengleich subgroups of a given\n");
      printf("space group up to conjugation in this space group!\n");
      printf("\n");
      printf("file1: Spacegroup R in standard affine form without translations\n");
      printf("file2: Pointgroup P of R with correct normalizer and formspace\n");
      printf("file3: GAP information about the conjugation classes of\n");
      printf("       maximal subgroups of the point group P of R.\n");
      printf("file4: (optional) presentation of P\n");
      printf("\n");
      printf("Options:\n");
      printf("-h : Give this help.\n");
      printf("-a : Calculate subgroups up to conjugation under\n");
      printf("     the affine normalizer of R.\n");
      printf("-f : Write representatives in files `file1`_i.\n");
      printf("-t : Omit the last dim R generators.\n");
      printf("-g : Do not use the database, use your own GAP-words.\n");
      printf("     (Only for experts!)\n");
      printf("-c : Output for constructing the database. Do not use!\n");
      printf("\n");
      printf("Cf.: TSupergroups\n\n");
      exit(1);
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

   if (is_option('g')){
      P = get_bravais(FILENAMES[1]);
      gapwords = mget_mat(FILENAMES[2], &no);

      if (FILEANZ > 3){
         presentation = mget_mat(FILENAMES[3], &panz);
         if (panz > 1){
            fprintf(stderr, "you should only give a single matrix as presention\n");
            exit(3);
         }
      }
      else{
         base = get_base(P);
         strong = strong_generators(base, P, TRUE);
         presentation = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
         OPT[0] = 0;
         presentation[0] = pres(strong, P, OPT);
         for (i = 0; i < P->dim; i++){
            free_mat(base[i]);
            free_bahn(strong[i]);
            free(strong[i]);
         }
         free(strong);
         free(base);
      }

      /* Untergruppen berechnen und ausgeben */
      /* =================================== */
      no_neu = no;
      if (no > 0){
         subgroups = tsubgroup(R, P, presentation[0], gapwords, &no_neu,
                               is_option('a'), is_option('c'));
         if (!is_option('c')){
            for (i = 0; i < no_neu; i++){
               sprintf(comment, "orbitlength: %d, pointgroup order: %d",
                       subgroups[i]->orbitlength, subgroups[i]->pointgrouporder);
               sprintf(file, "%s_%d", FILENAMES[0], i + 1);
               if (is_option('f'))
                  put_bravais(subgroups[i]->R, file, comment);
               else
                  put_bravais(subgroups[i]->R, 0, comment);
            }
         }
      }
   }
   else{
      subgroups = tsubgroup_db(R, is_option('a'), &no_neu);
      for (i = 0; i < no_neu; i++){
         sprintf(comment, "orbitlength: %d, pointgroup order: %d",
                 subgroups[i]->orbitlength, subgroups[i]->pointgrouporder);
         sprintf(file, "%s_%d", FILENAMES[0], i + 1);
         if (is_option('f'))
            put_bravais(subgroups[i]->R, file, comment);
            else
               put_bravais(subgroups[i]->R, 0, comment);
      }
   }


   /* Speicherfreigabe */
   /* ================ */
   if (is_option('t')){
      R->gen_no += (R->dim - 1);
   }
   free_bravais(R);
   if (is_option('g')){
      free_bravais(P);
      free_mat(presentation[0]);
      free(presentation);
      for (i = 0; i < no; i++){
         free_mat(gapwords[i]);
      }
      free(gapwords);
   }
   for (i = 0; i < no_neu; i++){
      free_TSubgroup_TYP(subgroups[i]);
   }
   if (no > 0)
      free(subgroups);

   /* Debugging */
   /* ========= */
   if (INFO_LEVEL & 12){
      fprintf(stderr,"write pointer_statistics\n");
      pointer_statistics(0,0);
   }

   exit(0);
}



