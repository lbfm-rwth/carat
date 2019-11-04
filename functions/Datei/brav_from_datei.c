#include "typedef.h"
#include "getput.h"
#include "datei.h"
#include "idem.h"

/***************************************************************************
@
@---------------------------------------------------------------------------
@
@ FILE: brav_from_datei.c
@
@---------------------------------------------------------------------------
@
****************************************************************************/


/***************************************************************************
@
@---------------------------------------------------------------------------
@
@ bravais_TYP *brav_from_datei(char *symb,int almost,int zclass)
@ 
@ Reads a single bravais group from the catalog.
@ The symbol, the number almost and zclass are exactly as in the
@ standalone Datei resp. Bravais_catalog
@
@---------------------------------------------------------------------------
@
****************************************************************************/
bravais_TYP *brav_from_datei(const char *symb,int almost,int zclass)
{

   bravais_TYP *RES;

   symbol_out *S;

   char *file;

   int i=1;

   /* initialize */
   S = read_symbol_from_string(symb);
   get_zentr(S);

   while (i<almost){
      /* get the next almost decomposable group */
      i++;
      file = S->fn;
      free_bravais(S->grp);
      free(S);
      S = get_symbol(file);
      if (file != NULL) free(file);
   }


   if (zclass == 1){
      RES = S->grp;
   }
   else{
      RES = Z_class(S->grp,S->grp->zentr[zclass-2]);
      free_bravais(S->grp);
   }

   if (S->fn != NULL) free(S->fn);
   free(S);

   return RES;

}

