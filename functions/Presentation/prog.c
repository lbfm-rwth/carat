
/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: prog.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

#include "typedef.h"
#include"matrix.h"
#include"symm.h"
#include"getput.h"
#include"bravais.h"
#include"sort.h"
#include"polyeder.h"
#include"presentation.h"
#include"tools.h"
#include"tietzetrans.h"

        extern char **FILENAMES;
        extern int FILEANZ;
	
int SFLAG;

static int lookfor(int a, int *b,int anz);
static int lookfor(int a, int *b,int anz);

static int lookfor(int a, int *b,int anz)
{
  int 		i;

  if(a > 0){
     for(i=0; i<anz; i++){
	if(a == b[i])
	   return(i+1);
     }
  }
  else
     for(i=0; i<anz; i++){
	if(-a == b[i])
	   return((-1)*(i+1));
     }
 return(-1);
}/** lookfor() **/

/**************************************************************************\
@--------------------------------------------------------------------------
@ presentation_TYP *prog(x,G,Form)
@
@ matrix_TYP *x: 	affine vector, starting point of algorithm ie.(0,0,0,1) 
@ bravais_TYP *G:       given space group with translationlattice Z^n
@ matrix_TYP *Form:     positive definite G-invariant form
@
@ claculates a presentation of the spacegroup G where the translations of
@ G have to be Z^(nx1)
@
@--------------------------------------------------------------------------
\**************************************************************************/
presentation_TYP *prog(x,G,Form)
matrix_TYP *x;
bravais_TYP *G;
matrix_TYP *Form;
{
  int			i,j;
  int			a;
  int 			anz;
  int 			**tmp1;
  int			s;
  int			ST_anz;
  int 			*Liste,*Liste_inv;
  int			*erzlist;
  matrix_TYP		*vec;
  polyeder_TYP 		*Pol;
  presentation_TYP 	*erg;
  anne_presentation_TYP *tmp;
  int			new_anz;
  int			**old_as_new;
  int			**new_as_old;

        extern polyeder_TYP *fub();
	extern anne_presentation_TYP *pres();
	extern void tilman_tietze();
  	extern int **back();

	Pol = fub(x,G,Form);

	vec = barzen(Pol);

        tmp = pres(Pol);
	ST_anz = tmp->gen_no;

	erg = make_pres(tmp);

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
   	   for(j=1; j <= tmp1[i][0]+1; j++){
       	      new_as_old[i][j] = lookfor(tmp1[i][j],erzlist, ST_anz);
	   }
 	   free(tmp1[i]);
    	}
	free(tmp1);

	old_as_new = (int**)malloc(ST_anz * sizeof(int*));
        for(i=0; i<ST_anz; i++){
	   j = erzlist[i];
	   old_as_new[i]= make_prod(Pol->wall[j-1]->word); 
	}

tilman_tietze(erg, new_anz, old_as_new, new_as_old,G->gen);

free_polyeder(Pol);

return(erg);
}
