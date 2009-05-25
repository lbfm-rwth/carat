#include "typedef.h"
#include "matrix.h"
#include "symm.h"
#include "getput.h"
#include "bravais.h"
#include "sort.h"
#include "polyeder.h"
#include "presentation.h"
#include "tools.h"
#include "tietzetrans.h"
#include "datei.h"

        extern char **FILENAMES;
        extern int FILEANZ;
	
int SFLAG;
int main (int argc, char *argv[])
{
  int			i,j;
  int			a;
  int			s;
  int 			anz;
  int			ST_anz;
  int 			**tmp1;
  int 			*Liste,*Liste_inv;
  int			*erzlist;
  matrix_TYP 		*x;
  matrix_TYP 		*Form;
  matrix_TYP		*vec;
  bravais_TYP 		*G;
  polyeder_TYP 		*Pol;
  presentation_TYP 	*erg;
  anne_presentation_TYP *tmp;
  char filename_out[100];
  int			new_anz;
  int			**old_as_new;
  int			**new_as_old;

        extern polyeder_TYP *fub();
	extern anne_presentation_TYP *pres();
	extern void tilman_tietze();
  	extern int **back();

  read_header(argc, argv);
  if(FILEANZ != 3)
  {
   printf("\n");
   printf("Usage\n");
   printf("  %s file1 file2 file3 \n", argv[0]);
   printf(" \n");
   printf(" where file1 contains a bravais_TYP describing\n");
   printf(" a set of affine matrices generating a space group.\n");
   printf(" \n");
   printf(" where file2 contains a matrix_TYP describing an affine vector \n");
   printf(" which is taken as a starting point for the algorithm, i.e. (0,0,0,1)\n");
   printf(" \n");
   printf(" file3 contains a matrix_TYP describing a positive definite form,\n");
   printf(" invariant under the point group of the given spacegroup\n");
   printf(" \n");
   printf(" Calculates a presentation of the group in file1, and writes\n");
   printf(" it to the file1.carat.\n");
   printf(" \n");
   printf(" The options are:\n");
   printf("\n");
   printf(" -G   : Output the presentation in Gap format on a file\n");
   printf("        called <file>.gap.short .\n");
   printf(" -h   : Gives you this help.\n");
   printf("\n");
   printf("\n");

   if (is_option('h')){
      exit(0);
   }
   else{
      exit(31);
   }
   }

        /* setting SFALG according to optionnumber('h') */
        if (is_option('h') && optionnumber('h') == 8){
           SFLAG = 1;
        }

        G = get_bravais(FILENAMES[0]);
	x = get_mat (FILENAMES[1]);
	Form = get_mat (FILENAMES[2]);

	Pol = fub(x,G,Form);

	vec = barzen(Pol);

        tmp = pres(Pol);
	ST_anz = tmp->gen_no;

	erg = make_pres(tmp);
/*
	printf("erg->rhsproduct \n");
	 for(i=0; i<erg->norelators; i++){
	     for(j=0; j<erg->relators[i].rhsnproduct; j++)
       		  printf(" %d",erg->relators[i].rhsproduct[j]);
	     printf("\n");
	 }
*/
	Liste = generate_Liste(Pol, &Liste_inv);
	erzlist = (int*)malloc(ST_anz * sizeof(int));
	s = 0;
	for(i=0; i<Pol->wall_no; i++){
	   if(Liste[i] > 0){
	      erzlist[s] = Liste[i];
	      s++;
	   }
 	}
	anz = 0;
 	tmp1 = back(vec,G->gen,G->gen_no,Pol,Form);
	/* umschreiben Seitentrafo's auf die Erzeuger (Liste[i]>0) */
	for(i=0; i<G->gen_no; i++){
	    for(j=1; j<=tmp1[i][0]; j++){ 
	  	a = tmp1[i][j];
		tmp1[i][j] = Liste[a-1];
	    }
	}
	
 	new_anz = G->gen_no;
	new_as_old = (int**)malloc(new_anz * sizeof(int*));
	for(i=0; i<G->gen_no; i++){
	   new_as_old[i] = (int*)malloc((tmp1[i][0]+1) * sizeof(int));
	   new_as_old[i][0] = tmp1[i][0];
   	   for(j=1; j <= tmp1[i][0]+1; j++)
       	      new_as_old[i][j] = tmp1[i][j];
	/*
   	   for(j=1; j <= tmp1[i][0]+1; j++){
 	      if(tmp1[i][j] > 0)
       	      new_as_old[i][j] = tmp1[i][j];
	      else
       	      new_as_old[i][j] = Liste_inv[(-1)*tmp1[i][j]]; 
    	   }
	*/
 	   free(tmp1[i]);
    	}
	free(tmp1);
/*
printf(" Woerter fuer G \n");
for(i=0; i< G->gen_no; i++){
for(j=0; j<new_as_old[i][0]+1; j++)
printf(" %d ",new_as_old[i][j]);
printf("\n ");
}
printf("\n ");
*/
	anz = 0;

	old_as_new = (int**)malloc(ST_anz * sizeof(int*));
        for(i=0; i<ST_anz; i++){
	   j = erzlist[i];
	   old_as_new[i]= make_prod(Pol->wall[j-1]->word); 
	}
/*
printf(" Relationen fuer ST \n");
for(i=0; i<ST_anz; i++){
for(j=0; j<old_as_new[i][0]+1; j++)
printf(" %d ",old_as_new[i][j]);
printf("\n ");
}
*/
/* tietzetrans ist nur fuer endliche Gruppen geschrieben. nehme noch die Inversen der Erzeuger hinzu, und aendere die Worte entsprechend ab. */

/* free_derivedsg(erg->generators);
   exit(2); */ 
/* -------------------------------------------------------------------------
    Transformieren einer gegebenen Praesentation presentation mit Erzeugern
    wie in grouggennew in eine Praesentation mit Erzeugern wie sie in
    groupgenold abgelegt sind. groupgenold beinhaltet ausserdem die Dar-
    stellung jedes Erzeugers als Produkt der Erzeuger in groupgennew.
    Ebenso beinhaltet groupgennew die Darstellung jedes Elementes als 
    Produkt der Elemente aus groupgenold. Diese Angaben werden wie unten
    erklaert fuer die Tietzetransformationen benoetigt.
----------------------------------------------------------------------------*/

tilman_tietze(erg, new_anz, old_as_new, new_as_old,G->gen);

fprintf(stderr,"\n Output files are: \n");

if (is_option('G')){
   sprintf(filename_out, "%s.gap.short", FILENAMES[0]);
   put_presentation(erg,filename_out,"G");
   printf("presnet nach tietze\n");
   put_presentation(erg,NULL,"G");
   put_presentation(erg,NULL,"C");
}

   sprintf(filename_out, "%s.carat", FILENAMES[0]);
   put_presentation(erg,filename_out,"C");

free_polyeder(Pol);
free_mat(x); x = NULL;
free_bravais(G);
free_mat(Form); Form = NULL;
        if (SFLAG == 1){
           pointer_statistics(0,0);
        }
   printf("\n");
}
