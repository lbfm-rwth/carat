#include "typedef.h"
/* #include "const.h" */
#include "dertype.h"
#include "tietzetrans.h"

extern derived_TYP* copy_derived();
extern derivedsg_TYP* copy_derivedsg();
extern void dumprelators();
extern int is_in_list();
extern int is_in_list_l();
extern void insert_list();
extern matrix_TYP* copy_mat();
extern matrix_TYP* mat_mul();
extern derivedsg_TYP* init_derivedcomplete();

/*-----------------------------------------------------------------------------
   Prozedur, die einen Relator aus der Praesentation loescht
------------------------------------------------------------------------------*/

void delete_relator(presentation, relatorindex)
presentation_TYP* presentation;
int relatorindex;
{
int i,j;
for (i=relatorindex; i<(presentation->norelators-1);i++)
    {
    if (presentation->relators[i].lhsnproduct>0)
       {
       free(presentation->relators[i].lhsproduct);
       presentation->relators[i].lhsnproduct=0;
       presentation->relators[i].lhsproduct=NULL;
       }
    presentation->relators[i].lhsnproduct=presentation->relators[i+1].lhsnproduct;
    if (presentation->relators[i].lhsnproduct!=0)
       {
       presentation->relators[i].lhsproduct = (int*) malloc (sizeof(int)*presentation->relators[i].lhsnproduct);
       for (j=0;j<presentation->relators[i].lhsnproduct;j++)
           presentation->relators[i].lhsproduct[j]=presentation->relators[i+1].lhsproduct[j];
       }
    if (presentation->relators[i].rhsnproduct>0)
        {
        free(presentation->relators[i].rhsproduct);
        presentation->relators[i].rhsnproduct=0;
        presentation->relators[i].rhsproduct=NULL;
        }
    presentation->relators[i].rhsnproduct=presentation->relators[i+1].rhsnproduct;
    if (presentation->relators[i].rhsnproduct!=0)
       {
       presentation->relators[i].rhsproduct = (int*) malloc (sizeof(int)*presentation->relators[i].rhsnproduct);
       for (j=0;j<presentation->relators[i].rhsnproduct;j++)
           presentation->relators[i].rhsproduct[j]=presentation->relators[i+1].rhsproduct[j];
       }
    }
presentation->norelators--;
if (presentation->relators[presentation->norelators].lhsnproduct>0)
   {
   free(presentation->relators[presentation->norelators].lhsproduct);
   presentation->relators[presentation->norelators].lhsproduct = NULL;
   presentation->relators[presentation->norelators].lhsnproduct = 0;
   }
if (presentation->relators[presentation->norelators].rhsnproduct>0)
   {
   free(presentation->relators[presentation->norelators].rhsproduct);
   presentation->relators[presentation->norelators].rhsproduct = NULL;
   presentation->relators[presentation->norelators].rhsnproduct = 0;
   }
}


/*-----------------------------------------------------------------------------
   Funktion, die den Relator sucht, der den Generator mit Index gen als
   Produkt von anderen Erzeugern darstellt.
------------------------------------------------------------------------------*/
int which_relator(presentation, gen)
presentation_TYP* presentation;
int gen;
{
int i;
for (i=0;i<presentation->norelators;i++)
    {
    if ((presentation->relators[i].lhsnproduct==1) && (presentation->relators[i].lhsproduct[0]==gen))
        {
        return(i);  
        }
    }
}

/*---------------------------------------------------------------------------
   Prozedur, die ein Produkt insproduct in ein anderes Produkt product an
   der Stelle position einfuegt. Das Element an position wird dabei 
   geloescht.
----------------------------------------------------------------------------*/
void insert_in_product(product, position, nproduct, insproduct, ninsproduct)
int** product;
int position;
int *nproduct;
int*insproduct;
int ninsproduct;
{
int naltproduct;
int i;
naltproduct=(*nproduct);
*nproduct = (*nproduct)+ninsproduct-1;
*product=(int*) realloc (*product, sizeof(int)*(*nproduct));
/* Elemente ab Position (position+1) werden an das Ende des nun groesseren
   Feldes geschoben */
for (i=naltproduct-2;i>=position;i--)
    {
    (*product)[i+ninsproduct]=(*product)[i+1];
    }
/* Die neuen Elemente werden an den freigewordenen Plaetzen eingefuegt */
for (i=0;i<ninsproduct;i++)
    {
    (*product)[i+position]=insproduct[i];
    }
}

/*---------------------------------------------------------------------------
   Prozedur, die eine bestimmte Anzahl von Eintraegen (number) in einem
   Integer-Feld loescht, ab einer bestimmten Position.
---------------------------------------------------------------------------*/
void remove_product(product, nproduct, position, number)
int** product;
int* nproduct;
int position;
int number;
{
int i;
for (i=position;i<(*nproduct);i++)
    {
    (*product)[i]=(*product)[i+number];
    }
*nproduct = (*nproduct)-number;
*product=(int*) realloc (*product, sizeof(int)*(*nproduct));
}

/*----------------------------------------------------------------------------
   Diese Funktion stellt fest, wie oft der Wert an der Position position
   hintereinander im Feld auftritt.
----------------------------------------------------------------------------*/
int block_detect(product, nproduct, position)
int * product;
int nproduct;
int position;
{
int blockcount=1;
while (((position+blockcount)<nproduct) && (product[position]==product[position+blockcount]))
     blockcount++;
return(blockcount);
}


/*-------------------------------------------------------------------------
   Prozedur, die ein Produkt an ein anderes Produkt anhaengt.
-------------------------------------------------------------------------*/
void add_product(product, nproduct, insproduct, ninsproduct)
int** product;
int *nproduct;
int*insproduct;
int ninsproduct;
{
int naltproduct;
int i;
naltproduct=(*nproduct);
*nproduct = (*nproduct)+ninsproduct;
*product=(int*) realloc (*product, sizeof(int)*(*nproduct));
for (i=0;i<ninsproduct;i++)
    {
    (*product)[i+naltproduct]=insproduct[i];
    }
}


/*---------------------------------------------------------------------------
   Diese Prozedur ersetzt in der gesamten Praesentation das auftauchen des
   Generators gen durch die rechte Seite des Relators relatorind. Dieser
   Relator sollte ein Darstellung des Generators als Produkt von anderen
   Generatoren beinhalten. Dies kann sichergestellt werden, durch vor-
   heriges aufrufen von which_relator.
---------------------------------------------------------------------------*/
void substitute_relator (presentation, relatorind, gen)
presentation_TYP* presentation;
int relatorind;
int gen;
{
int i, j;
for (i=0;i<presentation->norelators;i++)
    {
    j=0;
    while (j<presentation->relators[i].lhsnproduct)
        {
        if (presentation->relators[i].lhsproduct[j]==gen)
           {
           insert_in_product(&presentation->relators[i].lhsproduct, j, &presentation->relators[i].lhsnproduct, presentation->relators[relatorind].rhsproduct, presentation->relators[relatorind].rhsnproduct);
           j+=presentation->relators[relatorind].rhsnproduct;
           }
        else
           {
            j++;
           }
        }
    j=0;
    while (j<presentation->relators[i].rhsnproduct)
        {
        if (presentation->relators[i].rhsproduct[j]==gen)
           {
           insert_in_product(&presentation->relators[i].rhsproduct, j, &presentation->relators[i].rhsnproduct, presentation->relators[relatorind].rhsproduct, presentation->relators[relatorind].rhsnproduct);
           j+=presentation->relators[relatorind].rhsnproduct;
           }
        else
           {
           j++;
           }
        }
    }
}

void act_relators(presentation, oldindex,newindex)
presentation_TYP* presentation;
int oldindex;
int newindex;
{
int i, j;
for (i=0;i<presentation->norelators;i++)
    {
    for (j=0;j<presentation->relators[i].lhsnproduct;j++)
        {
        if (presentation->relators[i].lhsproduct[j]==oldindex)
           {
           presentation->relators[i].lhsproduct[j]=newindex;
           }
        }
    for (j=0;j<presentation->relators[i].rhsnproduct;j++)
        {
        if (presentation->relators[i].rhsproduct[j]==oldindex)
           {
           presentation->relators[i].rhsproduct[j]=newindex;
           }
        }
    }
}
/* -------------------------------------------------------------------------
    Transformieren einer gegebenen Praesentation presentation mit Erzeugern
    wie in grouggennew in eine Praesentation mit Erzeugern wie sie in
    groupgenold abgelegt sind. groupgenold beinhaltet ausserdem die Dar-
    stellung jedes Erzeugers als Produkt der Erzeuger in groupgennew.
    Ebenso beinhaltet groupgennew die Darstellung jedes Elementes als 
    Produkt der Elemente aus groupgenold. Diese Angaben werden wie unten
    erklaert fuer die Tietzetransformationen benoetigt.
----------------------------------------------------------------------------*/
void tietzetrans(groupgenold, groupgennew, presentation)
derivedsg_TYP *groupgenold, *groupgennew;
presentation_TYP* presentation;
{
int i,j,k;
int relatorind;
int index;
int oldindex;
int lauf;
int number;
int blockcount;
int** inverse_product;
int * inverse_nproduct;
int orderno;
matrix_TYP* hilf2;
derivedsg_TYP* hilfcomplete;
derived_TYP* hilf;
derived_TYP* hilfgen;
char filename[100];

/* Alle urspruenglichen Erzeuger der Gruppe in die Praesentation aufnehemen,
   falls sie noch nicht vorhanden sind. Gleichzeitig wird ein Relator hin-
   zugenommen, der das neue Element in der Generatorliste als Product der
   alten Elemente ausdrueckt (Tietzetranformation T3 (s. "Combinatorial
   group theory", Magnus-Karrass-Solitar, Dover Publications, Inc. 1976
   Seite 49f. ) */
for (i=0;i<groupgenold->firstfree;i++)
    {
    if (!is_in_list_l(presentation->generators, groupgenold->list[i],&index))
       {
       if (presentation->norelators>=(presentation->ext_factor*EXT_SIZE))
          {
          presentation->ext_factor++;
          presentation->relators=(relator_TYP*) realloc(presentation->relators, sizeof(relator_TYP)*EXT_SIZE*presentation->ext_factor);
          }
       hilfgen=copy_derived(groupgenold->list[i]);
       insert_list(presentation->generators, hilfgen);
      /*
     printf("\n im tietze ... \n");
      dumplist(presentation->generators);
      */
       k=presentation->norelators;
    /*  printf(" add generator %d + add relator %d  \n",presentation->generators->firstfree,presentation->norelators);*/
       presentation->relators[k].lhsnproduct = 1;
       presentation->relators[k].lhsproduct = (int*) malloc (sizeof(int));
       presentation->relators[k].lhsproduct[0] = presentation->generators->firstfree-1;
       presentation->relators[k].rhsnproduct = groupgenold->list[i]->nproduct;
       presentation->relators[k].rhsproduct = (int*) malloc (sizeof(int)*presentation->relators[k].rhsnproduct);
       for (j=0;j<presentation->relators[k].rhsnproduct;j++)
           {
           presentation->relators[k].rhsproduct[j] = groupgenold->list[i]->product[j];
            
           }
       presentation->norelators++;
       }
    }

if (is_option('d')){
   gap_out(presentation,"gg.1");
}
   
/* An diesem Punkt liegt eine Praesentation vor, die alle urspruenglichen 
   Generatoren beinhaltet (die gewollten), aber zusaetzlich die bei der
   Kommutatorreihenbestimmung generierten. Diese muessen nun eleminiert 
   werden. Dazu werden zunaechst neue Relatoren hinzugefuegt, und zwar 
   diejenigen, die die "ungewollten" Generatoren als Produkt der "gewollten"
   darstellen. Diese liegen in anderer Form in groupgennew bereits vor.
   Hierbei handelt es sich um ein Tietzetransformation T1 (s.o.)*/
/*
printf(" Status vor T1: \n");
dumprelators(presentation->relators,presentation->norelators);
*/
for (i=0;i<groupgennew->firstfree;i++)
    {
    if (!is_in_list_l(groupgenold, groupgennew->list[i],&index))
       {
 /*     printf(" add relator %d  \n",presentation->norelators);*/
       k=presentation->norelators;
       presentation->relators[k].lhsnproduct = 1;
       presentation->relators[k].lhsproduct = (int*) malloc (sizeof(int));
       presentation->relators[k].lhsproduct[0] = i;
       presentation->relators[k].rhsnproduct = groupgennew->list[i]->nproduct;
       presentation->relators[k].rhsproduct = (int*) malloc (sizeof(int)*(presentation->relators[k].rhsnproduct +1));
       for (j=0;j<presentation->relators[k].rhsnproduct;j++)
           {
           is_in_list_l(presentation->generators, groupgenold->list[groupgennew->list[i]->product[j]], &index);
           presentation->relators[k].rhsproduct[j] = index;
           }
       presentation->norelators++;
       } 
    }

if (is_option('d')){
   gap_out(presentation,"gg.1");
}
   
/* Nun ist die Praesentation in einem Status, in dem alle ungewollten 
   Generatoren einen entsprecheneden Relator haben, der sie als Produkt der
   gewollten ausdrueckt. Damit sind die Vorraussetzungen fuer die Tietze-
   transformation T4 erfuellt und alle ungewollten Erzeuger koennen 
   geloescht werden, ebenso die Relatoren, die sie als Produkt der anderen
   Erzeuger ausdruecken. Allerdings muessen auch in allen anderen Relatoren
   fuer die geloeschten Elemente das Produkt aus den verbleibenden Erzeugern 
   substituiert werden. */

for (i=0;i<groupgennew->firstfree;i++)
    {
    if (!is_in_list_l(groupgenold, groupgennew->list[i],&index))
       {
       is_in_list_l(presentation->generators, groupgennew->list[i], &index);
       relatorind=which_relator(presentation, index);
       substitute_relator (presentation, relatorind, index); 
       delete_relator(presentation, relatorind);
       }
    }
/* Ersetze die Indizes in den Relatoren (bisher bezueglich generators) durch Indizes bezueglich oldgen. Erspart umstaendliches loeschen in der Generatorliste. */

if (is_option('d')){
   gap_out(presentation,"gg.1");
}
   
for (i=0;i<presentation->generators->firstfree;i++)
   {
   act_relators(presentation,i,i+presentation->generators->firstfree);
   }
for (i=0;i<groupgenold->firstfree;i++)
    {
    is_in_list_l(presentation->generators,groupgenold->list[i],&oldindex);
    act_relators(presentation,oldindex+presentation->generators->firstfree,i);
    }
free_derivedsg(presentation->generators);
presentation->generators=copy_derivedsg(groupgenold);
/*
cayley_out(presentation, "group.t42");
*/
/*printf(" Status vor inverses... T4: \n");
dumprelators(presentation->relators,presentation->norelators);*/
/* Nun liegt eigentlich die gewuenschte Praesentation vor. Sie muss aller-
   dings noch auf die richtige Form gebracht werden, d.h. alle Relatoren
   muessen die Form haben: R(a...z)=1. Zur Zeit haben wir R(a..z)=R(a..z).
   Es muss als jeder Relator solange von rechts mit Inversen der Generatoren
   multipliziert werden, bis auf der rechten Seite nur noch die Eins uebrig
   bleibt.
*/
/*
cayley_out(presentation, "group.t43");
*/
/* Bestimme Wort fuer das Inverse eines jeden Generators */
inverse_product=(int**) malloc(sizeof(int*)*presentation->generators->firstfree);
inverse_nproduct=(int*) malloc(sizeof(int)*presentation->generators->firstfree);
hilfcomplete=init_derivedcomplete(presentation->generators->list[0]->element->cols);
hilf=(derived_TYP*) malloc (sizeof(derived_TYP));
hilf->nproduct = 0;
for (i=0;i<presentation->generators->firstfree;i++)
    {
    orderno=1;
    hilf->element=copy_mat(presentation->generators->list[i]->element);
    do
       {
       orderno++;
       hilf2=copy_mat(hilf->element);
       free_mat(hilf->element);
       hilf->element=mat_mul(hilf2,presentation->generators->list[i]->element);
       free_mat(hilf2);
       } while (!is_in_list(hilfcomplete,hilf,&index));    
    inverse_product[i]=(int*) malloc (sizeof(int)*(orderno-1));
    inverse_nproduct[i]=orderno-1;
    for (j=0;j<orderno-1;j++)
        {
        inverse_product[i][j]=i;
        }
    }
free_derived(hilf);
free_derivedsgcomp(hilfcomplete);
for (i=0;i<presentation->norelators;i++)
    {
    while (presentation->relators[i].rhsnproduct>0)
         {
         j=presentation->relators[i].rhsproduct[presentation->relators[i].rhsnproduct-1];
         add_product(&presentation->relators[i].lhsproduct,&presentation->relators[i].lhsnproduct, inverse_product[j], inverse_nproduct[j]);          
         presentation->relators[i].rhsnproduct--;
         }
    if (presentation->relators[i].rhsnproduct>0) 
       {
       free(presentation->relators[i].rhsproduct);
       presentation->relators[i].rhsnproduct=0;
       presentation->relators[i].rhsproduct=NULL;
       }
    }
/*
cayley_out(presentation, "group.t5");
*/
/* Die Praesentation ist nun in der richtigen Form usw.. Sie muss jetzt nur
   noch gekuerzt werden, falls es moeglich ist. Dazu:
     1. Falls ein Element n mal hintereinander auftaucht im Produkt und die
        Ordnung m hat, wird das Element nur n mod m mal aufgefuehrt.
     2. Falls ein Relator leer ist wird er geloescht
*/
for (i=0;i<presentation->norelators;i++)
    {
    lauf=0;
    while (lauf<presentation->relators[i].lhsnproduct)
       {
       blockcount=block_detect(presentation->relators[i].lhsproduct,presentation->relators[i].lhsnproduct, lauf);
       if ((blockcount>
          inverse_nproduct[presentation->relators[i].lhsproduct[lauf]])
        &&(!((blockcount==presentation->relators[i].lhsnproduct)
           &&(blockcount
          ==inverse_nproduct[presentation->relators[i].lhsproduct[lauf]]+1))))
          {
          number=0;
          while (number<=blockcount)
                number+=(inverse_nproduct[presentation->relators[i].lhsproduct[lauf]]+1);
          number-=(inverse_nproduct[presentation->relators[i].lhsproduct[lauf]]+1);
          remove_product(&presentation->relators[i].lhsproduct,
                &presentation->relators[i].lhsnproduct,lauf,number);
          lauf=0;
          }
       else
          lauf+=blockcount;
       }
    }
i=0;
while (i<presentation->norelators)
    {
    if (presentation->relators[i].lhsnproduct==0)
       delete_relator(presentation, i);
    else
       i++;
    }
for (i=0;i<presentation->generators->firstfree;i++)
    {
    free(inverse_product[i]);
    }
free(inverse_product);
free(inverse_nproduct);
}
