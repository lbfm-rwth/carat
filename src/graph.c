/* author: Oliver Heidbuechel */
/* last change: 16.03.2001 */


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



int INFO_LEVEL;
extern int SFLAG;
boolean GRAPH_DEBUG = FALSE;
boolean GRAPH = FALSE;

int main (int argc, char *argv[])
{
   matrix_TYP **presentation,
               *erg,
              **base;

   bahn **strong;

   bravais_TYP *G;

   Q_data_TYP *data;

   int i,
       panz,
       OPT[6];



   read_header(argc, argv);
   if (FILEANZ < 1 || FILEANZ > 2 || (is_option('h') && optionnumber('h') == 0)){
      printf("\n");
      printf("Usage: %s 'file1' ['file2'] -[options]\n",argv[0]);
      printf("\n");
      printf("file1: REDUCED pointgroup G with CORRECT order\n");
      printf("file2: (Optional) Presentation of G\n");
      printf("\n");
      printf("Calculates the \"graph of inclusions\" for a\n");
      printf("geometric class given by G.\n");
      printf("For further information on the output see\n");
      printf("$CARATPATH/tex/Graph\n");
      printf("\n");
      printf("Options:\n");
      printf("-i    : Print the Z-classes to 'file1.i' and the affine classes\n");
      printf("        to 'file1.i.j',\n");
      printf("-f    : Calculates the formspace even if it is given.\n");
      printf("        CAUTION:\n");
      printf("        If you give the formspace and the normalizers,\n");
      printf("        they have to be correct. Then G->gen and G->normal\n");
      printf("        have to generate the normalizer of G.\n");
      printf("-l    : If the order of a cohomology group is >= %i,\n", TWOTO21);
      printf("        you have to use this option; the program needs more\n");
      printf("        time with this option.\n");
      printf("-d    : only for debugging\n");
      exit(11);
   }

   INFO_LEVEL = optionnumber('h');
   if (INFO_LEVEL & 12){
      SFLAG = 1;
   }

   /* get data */
   G = get_bravais(FILENAMES[0]);
   if (FILEANZ == 2){
      presentation = mget_mat(FILENAMES[1],&panz);
      if (panz > 1){
         fprintf(stderr, "you should only give a single matrix as presention\n");
         exit(12);
      }
   }
   else{
      base = get_base(G);
      strong = strong_generators(base,G,TRUE);
      presentation = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
      OPT[0] = 0;
      presentation[0] = pres(strong, G, OPT);
   }
   if (G->form == NULL || G->form_no == 0 || is_option('f')){
     if (G->form != NULL){
        for (i=0; i< G->form_no; i++)
           free_mat(G->form[i]);
        free(G->form);
     }
     G->form = formspace(G->gen,G->gen_no,1, &G->form_no);
   }
   if (is_option('d')){
      GRAPH_DEBUG = TRUE;
   }
   data = get_Q_data(G, presentation[0], is_option('l'));

   /* write informations about the Q-class */
   put_Q_data(data, FILENAMES[0], is_option('i'));

   /* calculate the graph for subgroups */
   erg = subgroupgraph(data);
   put_mat(erg,0,"inclusions for all spacegroups",0);


   /* clean up  */
   free_mat(erg);
   free_Q_data(data);
   free_mat(presentation[0]);
   free(presentation);
   if (FILEANZ == 1){
      for (i = 0; i < G->dim; i++){
         free_mat(base[i]);
         free_bahn(strong[i]);
         free(strong[i]);
      }
      free(strong);
      free(base);
   }
   free_bravais(G);

   /* for debugging */
   if (INFO_LEVEL & 12){
      fprintf(stderr,"write pointer_statistics\n");
      pointer_statistics(0,0);
   }

   exit(0);
}
