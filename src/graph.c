/* author: Oliver Heidbuechel */
/* last change: 14.09.2000 */


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
boolean GRAPH_DEBUG;


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
       OPT[6],
       U[6], O[6];



   read_header(argc, argv);
   if (FILEANZ < 1 || FILEANZ > 2 || (is_option('h') && optionnumber('h') == 0)){
      printf("\n");
      printf("Usage: %s 'file1' [-h] [-o] [-i] [-f] ['file2']\n",argv[0]);
      printf("\n");
      printf("file1: REDUCED pointgroup G with CORRECT order\n");
      printf("file2: (Optional) Presentation of G\n");
      printf("\n");
      printf("Calculates the graph of inclusions of the Q-class given by G.\n");
      printf("For further information on the output see example 13 of\n");
      printf("the CARAT introduction (http://wwwb.math.rwth-aachen.de/carat).\n");
      printf("\n");
      printf("Options:\n");
      printf("-h    : Give this help.\n");
      printf("-o    : Do not calculate the corresponding supergroup numbers.\n");
      printf("        The programm is faster then.\n");
      printf("-i    : Print the Z-classes to 'file1.i' and the affine classes\n");
      printf("        to 'file1.i.j',\n");
      printf("-f    : Recalculate the formspace even if it is given.\n");
      printf("-d    : Only for debugging!\n");
      printf("\n");
      printf("CAUTION: If the formspace and the normalizer are given,\n");
      printf("         they have to be correct.\n");
      printf("\n");
      printf("Cf.: KSupergroups, KSubgroups\n");
      printf("\n");
      exit(11);
   }

   INFO_LEVEL = optionnumber('h');
   if (INFO_LEVEL & 12){
      SFLAG = 1;
   }

   /* get data */
   G = get_bravais(FILENAMES[0]);


   /* trivial cases */
   if (G->order == 0){
      printf("There is 1 Z-Class with 1 Space Group!\n");
      erg = init_mat(1,1,"");
      put_mat(erg,0,0,0);
      free_mat(erg);
      free_bravais(G);
      exit(0);
   }
   if (G->order == 2){
      for (i = 0 ; i < G->gen_no; i++){
         if (!G->gen[i]->flags.Scalar){
	    break;
	 }
      }
      if (i == G->gen_no){
         U[0] = 2; U[1] = 6; U[2] = 14; U[3] = 30; U[4] = 62; U[5] = 126;
         O[0] = 1; O[1] = 3; O[2] = 7; O[3] = 15; O[4] = 31; O[5] = 63;
         printf("There is 1 Z-Class with 1 Space Group!\n");
         printf("1: 1 (%i", U[G->dim - 1]);
         if (!is_option('o'))
	    printf(", %i", O[G->dim - 1]);
	 printf(", 2^1)\n");
         erg = init_mat(1,1,"1");
         put_mat(erg,0,0,0);
         free_bravais(G);
         free_mat(erg);
         exit(0);
      }
   }

   /* get more data */
   if (FILEANZ == 2){
      presentation = mget_mat(FILENAMES[1],&panz);
      if (panz > 1){
         fprintf(stderr, "You should only give a single matrix as presention!\n");
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

   /* calculate the graph */
   erg = subgroupgraph(data, !is_option('o'));
   put_mat(erg,0,0,0);


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






