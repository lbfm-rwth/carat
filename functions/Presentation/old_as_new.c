#include "typedef.h"
#include "dertype.h"
#include "tietzetrans.h"

extern matrix_TYP* mat_mul();
extern matrix_TYP* init_mat();
extern matrix_TYP* copy_mat();
extern void free_mat();
extern matrix_TYP* mat_inv();
/* Funktionen aus derived.c, zur Verwaltung der Element und Erzeugerlisten*/
extern int is_in_list();
extern void insert_list();
extern derivedsg_TYP* init_derivedcomplete();
extern derivedsg_TYP* init_derivedgen();
extern derivedsg_TYP* copy_derivedsg();
extern int* copy_product();
extern void free_derivedsg();
extern void free_derivedsgcomp();




void old_as_new(groupgenold, groupgennew, dim)
derivedsg_TYP* groupgenold;
derivedsg_TYP* groupgennew;
int dim;
{
derived_TYP* hilf;
int i,j,k;
int l;
int start, anzahl;
int index;
char c;
int gencount;
derivedsg_TYP* groupcomplete;
/* printf(" SCOUNT: %d\n", SCOUNT);*/
groupcomplete = init_derivedcomplete(dim);
anzahl=groupcomplete->firstfree;
start=0;
/* Die neuen Elemente der Liste befinden sich immer zwischen start und
   anzahl.
 */
gencount = 0;
while ((start!=anzahl)&&(gencount<groupgennew->firstfree)) 
   {
   for (i=start;i<anzahl;i++)
       {
       for (k=0;k<groupgennew->firstfree;k++)
           {
           hilf=(derived_TYP *)malloc(sizeof(derived_TYP));
	   hilf->nproduct = 0;
           hilf->element=mat_mul(groupcomplete->list[i]->element,groupgennew->list[k]->element); 
           if (!is_in_list(groupcomplete,hilf,&index))
             {
             hilf->nproduct=groupcomplete->list[i]->nproduct+1;
             hilf->product=(int*)malloc(sizeof(int)*
                                hilf->nproduct);
             for (j=0;j<groupcomplete->list[i]->nproduct;j++)
                 hilf->product[j]=groupcomplete->list[i]->product[j];
             hilf->product[hilf->nproduct-1]=k;
             insert_list(groupcomplete,hilf);   
             for (l=0;l<groupgenold->firstfree;l++)
                 {
                 if (is_in_list(groupcomplete, groupgenold->list[l],&index))
                    {
                    if (groupgenold->list[l]->nproduct == 0)
                       {
                       groupgenold->list[l]->product=copy_product(groupcomplete->list[index]);
                       groupgenold->list[l]->nproduct=groupcomplete->list[index]->nproduct;
                       gencount++;
                       }
                    }
                 }
             }
          else
             {
             free_derived(hilf);
             }
          }
       } 
   start=anzahl;
   anzahl=groupcomplete->firstfree;
   }
/* printf(" SCOUNT: %d\n", SCOUNT);*/
free_derivedsgcomp(groupcomplete);
/* printf(" SCOUNT: %d\n", SCOUNT);*/
}
