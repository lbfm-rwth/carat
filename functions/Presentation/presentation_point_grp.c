
/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: presentation_point_grp.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

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
#include "presentation.h"

/**************************************************************************\
@--------------------------------------------------------------------------
@ presentation_TYP *presentation_point_grp(bravais_TYP *G)
@
@ calculates a presentation of a finite integral matrix group G
@
@--------------------------------------------------------------------------
\**************************************************************************/
presentation_TYP *presentation_point_grp(bravais_TYP *G)
{
  int		i,j,s;
  int		c,r;
  int		n;
  matrix_TYP	*I;
  matrix_TYP	*F;
  matrix_TYP	*Y;
  matrix_TYP	**ttmp;
  bravais_TYP 	*H;
  presentation_TYP *erg;
  int		new_anz;
  int		**old_as_new;
  int		**new_as_old;
  char 		filename_out[100];
  
	extern presentation_TYP *prog();
	extern void tilman_tietze();

  n = G->dim;
  I = einheitsmatrix(G->dim);
  Y = init_mat(n+1,1,"");
  Y->array.SZ[n][0] = 1;

  H = copy_bravais(G);

 	c = G->gen[0]->cols+1;
	r = G->gen[0]->rows+1;

           j = G->gen_no;
	   ttmp = (matrix_TYP **)calloc(G->gen_no, sizeof(matrix_TYP *));
	   for (i=0; i<j ;i++){
	        ttmp[i] = copy_mat(H->gen[i]);
	        Check_mat(ttmp[i]);
	   }
	   F = rform(ttmp,j,I,101);
           for (i=0; i<j ;i++){
	       free_mat(ttmp[i]); ttmp[i] = NULL;
           }
	   free(ttmp);

	H->gen_no = H->gen_no + G->dim;
	H->gen = (matrix_TYP **)realloc(H->gen, H->gen_no *sizeof(matrix_TYP*));

        /*umsortieren*/
        for(i=1; i<=j; i++){
	    H->gen[n+j-i] = H->gen[j-i];
	    real_mat(H->gen[n+j-i],r,c);
	    H->gen[n+j-i]->array.SZ[r-1][c-1] = 1;
         } 

	/* only point_grp is relevant*/
        for (i=0; i<j ;i++){
	    for(s=0; s < r-1; s++)
	        H->gen[(G->dim)+i]->array.SZ[s][c-1] = 0;
           Check_mat(H->gen[(G->dim)+i]);
 	}	

	/*Translationen einfuegen*/
	for(i=0; i<n; i++){
	   H->gen[i] = einheitsmatrix(n+1);
	   H->gen[i]->array.SZ[i][c-1] = 1;
           Check_mat(H->gen[i]);
	}

	H->dim = H->gen[0]->rows;

	erg = prog(Y,H,F);

	/* omit the first n generators */
        old_as_new = (int**)malloc(H->gen_no * sizeof(int*));
        new_as_old = (int**)malloc((H->gen_no - n) * sizeof(int*));

	for(i=0; i<n; i++){
	    old_as_new[i] = (int*)malloc(sizeof(int));
	    old_as_new[i][0] = 0;
	}
	for(i=n; i<H->gen_no; i++){
	    old_as_new[i] = (int*)malloc(2*sizeof(int));
	    old_as_new[i][0] = 1;
	    old_as_new[i][1] = i-n+1;
	}
	for(i=0; i<H->gen_no-n; i++){
	   new_as_old[i] = (int*)malloc(2*sizeof(int));
	    new_as_old[i][0] = 1;
	    new_as_old[i][1] = n+i+1;
	}
 	new_anz = H->gen_no - n;
	memmove(H->gen,H->gen+n,new_anz * sizeof(matrix_TYP *));
	
	H->gen = (matrix_TYP **)realloc(H->gen, new_anz *sizeof(matrix_TYP*));
	H->gen_no = new_anz;
	
	tilman_tietze(erg,new_anz,old_as_new,new_as_old, H->gen);

if (is_option('G')){
   sprintf(filename_out, "%s.gap.short", FILENAMES[0]);
   put_presentation(erg,filename_out,"G");
   fprintf(stderr,"\n Output files are: file.gap.short \n");
}

/*
if (is_option('C')){
   sprintf(filename_out, "%s.carat", FILENAMES[0]);
   put_presentation(erg,filename_out,"C");
   fprintf(stderr,"\n Output files are: file.carat \n");
} */

free_bravais(H);
free_mat(Y); Y = NULL;
free_mat(F); F = NULL;

return erg;
}/*presentation_point_grp()*/
