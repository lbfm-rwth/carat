#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <longtools.h>
#include <presentation.h>
#include <sort.h>

/*************************************************************************
@
@-------------------------------------------------------------------------
@
@ void check_base(bahn **s,
@                 bravais_TYP *G)
@
@ checks the integrity of the structure s for the finite
@ group G, ie. checks all words in s[i]->words (0<=i<G->dim)
@ wheter they produce the right matrix if one insert G->gen
@ into them.
@ Neither s nor G are changed, and if an error occurs
@ the function gives a message to stderr.
@
@-------------------------------------------------------------------------
@
**************************************************************************/
void check_base(bahn **s,
                bravais_TYP *G){

   int i,
       j;


   matrix_TYP *M,
             **GENINV;

   GENINV = (matrix_TYP **) calloc(G->gen_no , sizeof(matrix_TYP *));

   for (i=0;i<G->dim;i++){
      for (j=0;j<s[i]->length;j++){
         if (s[i]->words[j] == NULL){
             fprintf(stderr,"error: no word\n");
         }
      }
   }


   for (i=0;i<G->dim;i++){
      for (j=0;j<s[i]->length;j++){
	 M = mapped_word(s[i]->words[j],G->gen,GENINV);
         if (mat_comp(M,s[i]->representatives[j])){
	   fprintf(stderr,"error: word is wrong i=%d j=%d\n",i,j);
	   put_word(s[i]->words[j],"G");
         }
         else if (FALSE){
	   fprintf(stderr,"word is good i=%d j=%d\n",i,j);
	   put_word(s[i]->words[j],"G");
         }
         free_mat(M);
      }
   }

   for (i=0;i<G->gen_no;i++)
     if (GENINV[i]) free_mat(GENINV[i]);

   if (GENINV) free(GENINV);

   return;

} /* check_base(.....) */

