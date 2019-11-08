#include "typedef.h"
#include "datei.h"

static int smaller(const char *s1,const char *s2)
{

   int i,
       j;

   /* the first criterion is the dimension */
   sscanf(s1,"%d",&i);
   sscanf(s2,"%d",&j);

   if  (i<j){
      return TRUE;
   }
   else if (i>j){
      return FALSE;
   }

   /* are these multiples of the same symbol ? */
   i = strcspn(s1,",");
   j = strcspn(s2,",");
   if ( (i==j) && strncmp(s1,s2,i)==0){
      if (strlen(s1)<strlen(s2)){
         return TRUE;
      }
      else{
         return FALSE;
      }
   }

   /* we are now in the position to distinguish it on one atomic symbol */
   /* the symbol has to contain a '-' now */
   s1 = strchr(s1,'-');
   s2 = strchr(s2,'-');
   s1++;
   s2++;
   sscanf(s1,"%d",&i);
   sscanf(s2,"%d",&j);

   if (i>j){
      return TRUE;
   }
   else if (i<j){
      return FALSE;
   }

   /* now one of them is of the form a-b', and the other one of a-b */
   /* a-b' < a-b */
   if (strchr(s1,'\'') == NULL){
      return TRUE;
   }

   return FALSE;
}

void right_order(char *string)
{

  int i,
      ordered,
      hom_no;

  char *irr_symbol[MAXDIM],
       *tmp,
       *tmp2,
       *tmp3;

  for (i=0;i<MAXDIM;i++){
     irr_symbol[i] = (char *) calloc(20 , sizeof(char));
  }

  tmp = (char *) calloc(20 * MAXDIM, sizeof(char));

  /* get the irreducible symbols, ie. those seperated by `;' */
  strcpy(tmp,string);
  tmp3 = tmp2 = tmp;
  hom_no = 0;
  while (tmp3 != NULL){
     tmp3 = strchr(tmp2,';');
     if (tmp3 != NULL)
        *tmp3 = 0; // cut the string at the semicolon
     strcpy(irr_symbol[hom_no],tmp2);
     if (tmp3 != NULL)
        tmp2 = tmp3+1; // if there was a semicolon, look at the text after it
     hom_no++;
  }


  /* order the symbols */
  ordered = FALSE;
  while (!ordered){
     ordered = TRUE;
     for (i=0;i<hom_no-1;i++){
        if (smaller(irr_symbol[i],irr_symbol[i+1])){
           tmp2 = irr_symbol[i];
           irr_symbol[i] = irr_symbol[i+1];
           irr_symbol[i+1] = tmp2;
           ordered = FALSE;
        }
     }
  }

  /* reprint them into string */
  sprintf(string,"%s",irr_symbol[0]);
  for (i=1;i<hom_no;i++){
     sprintf(tmp,"%s;%s",string,irr_symbol[i]);
     sprintf(string,"%s",tmp);
  }

  /* free */
  for (i=0;i<MAXDIM;i++){
     free(irr_symbol[i]);
  }
  free(tmp);

  return;

}
