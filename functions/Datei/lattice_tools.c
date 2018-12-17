#include "typedef.h"
#include "datei.h"
#include "getput.h"
#include "matrix.h"

/***************************************************************************
@
@---------------------------------------------------------------------------
@
@ FILE: lattice_tools.c
@
@---------------------------------------------------------------------------
@
***************************************************************************/


/***************************************************************************
@
@---------------------------------------------------------------------------
@
@ lattice_element *init_lattice_element()
@
@ initializes a lattice_element, and allocates memory for res->symbol
@ via calloc(100*sizeof(char));
@
@---------------------------------------------------------------------------
@
***************************************************************************/
lattice_element *init_lattice_element()
{
   lattice_element *res;

   res = (lattice_element *) malloc(1 * sizeof(lattice_element));

   res->grp = NULL;
   res->symbol = (char *) calloc(100 , sizeof(char));
   res->almost = -1;
   res->zclass = -1;
   res->alpha = -1;
   res->N_orbits = 0;
   res->TR = NULL;

   return res;
}

/***************************************************************************
@
@---------------------------------------------------------------------------
@
@ void free_lattice_element(lattice_element *x)
@
@ frees all stuff assigned in x, and x itself
@
@---------------------------------------------------------------------------
@
***************************************************************************/
void free_lattice_element(lattice_element *x)
{
   int i;

   if (x->grp != NULL) free_bravais(x->grp);
   for (i=0;i<x->N_orbits;i++) free_mat(x->TR[i]);
   free(x->symbol);
   if (x->TR != NULL) free(x->TR);

   free(x);
}

/***************************************************************************
@
@---------------------------------------------------------------------------
@
@---------------------------------------------------------------------------
@
***************************************************************************/
lattice_element *fget_lattice_element(FILE *F,int OPTION)
{
   lattice_element *E;
   int c;

   E = init_lattice_element();
   c=fscanf(F,"%s %d %d\n",E->symbol,&E->almost,&E->zclass);

   c=fscanf(F,"%d\n",&E->alpha);

   E->TR = fmget_mat(F,&E->N_orbits);

   /* sometimes we want to get hold of the group in the catalog as well */
   if (OPTION){
      E->grp = brav_from_datei(E->symbol,E->almost,E->zclass);
   }

   return E;
}

/***************************************************************************
@
@---------------------------------------------------------------------------
@
@---------------------------------------------------------------------------
@
***************************************************************************/
void fput_lattice_element(lattice_element *E,FILE *F)
{
   int i;

   if (F==NULL) F = stdout;

   fprintf(F,"Symbol: %s almost: %d zclass: %d\n",
              E->symbol,E->almost,E->zclass);

   fprintf(F,"alpha : %d\n",E->alpha);

   fprintf(F,"#%d\n",E->N_orbits);

   for (i=0;i<E->N_orbits;i++){
      fput_mat(F,E->TR[i],NULL,2);
   }

   if (E->grp != NULL) put_bravais(E->grp,NULL,NULL);

   return ;
}

