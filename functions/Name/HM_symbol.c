#include "typedef.h"
#include "name.h"
#include "tools.h"
#include "matrix.h"
#include "getput.h"
#include "bravais.h"
#include "contrib.h"
#include "sort.h"
#include "datei.h"



static void remove_underscore(char *s)
{

   char *PP;

   PP = strchr(s,'_');
   while (PP){
      *PP=' ';
      PP = strchr(PP,'_');
   }

}

/*************************************************************************
@ void display_HM_symbol(const char *qname,
@                        int zname1,
@                        int zname2,
@                        MP_INT *aff_name)
**************************************************************************/
void display_HM_symbol(const char *qname,
                       int zname1,
                       int zname2,
                       MP_INT *aff_name)
{

   int number,
       i,
       z1in,
       z2in,
       found=0;

   FILE *F;

   char qin[1028],
        affin[128],
        affstring[128],
        HMSYMBOL[128],
        hmtab[1024];;

   mpz_get_str(affstring,10,aff_name);

   sprintf(hmtab, "%s/qcatalog/translation_HM_symbol", get_data_dir());
   F = fopen(hmtab, "r");

   if (F == NULL){
      fprintf(stderr,"can't open %s\n", hmtab);
      exit(4);
   }

   fscanf(F,"#%d\n",&number);

   fprintf(stdout,"possible Herman-Mauguin symbols describing a group\n");
   fprintf(stdout,"isomorphic to the given one: \n");

   for (i=0;i<number;i++){
      fscanf(F,"qname: %s zname: %d %d aff_name: %s %s\n",
                      qin,&z1in,&z2in,affin,HMSYMBOL);

      if (strcmp(qin,qname) == 0
         && z1in == zname1
         && z2in == zname2
         && strcmp(affstring,affin) == 0){

         if (found > 0) fprintf(stdout," or ");
         remove_underscore(HMSYMBOL);
         fprintf(stdout,"%s",HMSYMBOL);
         found++;
      }

   }

   fprintf(stdout,"\n");

   if (found == 0){
      fprintf(stderr,"error in display_HM_symbol\n");
      exit(3);
   }

   return;
}

