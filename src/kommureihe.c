#include "typedef.h"
#include "dertype.h"
#include "tietzetrans.h"
#include "getput.h"
#include <stdio.h>
#define MAX_SERLEN 10

extern bravais_TYP* get_bravais();
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


/*------------------------------------------------------------------------ 
   Funktion, die ein Element mit einem anderen konjugiert und das Ergebnis
   zurueckliefert. Konjugiert wird nach ab(a^-1)
 ------------------------------------------------------------------------*/
derived_TYP* konjugate(konjelement, exponent)
derived_TYP* konjelement;
matrix_TYP* exponent;
{
matrix_TYP* hilf1, *hilf2;
derived_TYP* hilf;
hilf=(derived_TYP*)malloc(sizeof(derived_TYP));
hilf->nproduct = 0;
hilf->product = NULL;
hilf->left = NIL;
hilf->right = NIL;
hilf1=mat_mul(exponent, konjelement->element);
hilf2=mat_inv(exponent);
hilf->element=mat_mul(hilf1,hilf2);
free_mat(hilf1);
free_mat(hilf2);
return(hilf);
}


/*--------------------------------------------------------------------------
   Generate_sg, liefert das Erzeugnis, der in derivedgen abgelegten Gruppen-
   elemente. Das Erzeugnis wird dann in derivedcomplete abgelegt.
   Benutzte Funktionen:
      is_in_list  -- ueberprueft, ob eine bestimmtes Element in einer Liste
                     (beides Parameter) bereits enthalten ist. Der dritte
                     Parameter gibt die Stelle in der Liste an, an der es
                     gefunden wurde. Rueckgabewert ist 1 wenn das Element
                     in der Liste ist. Sonst 0.
     insert_list --  fuegt ein Element in eine Liste ein.
----------------------------------------------------------------------------*/
void generate_sg(derivedgen, derivedcomplete)
derivedsg_TYP* derivedgen, *derivedcomplete;
{
derived_TYP* hilf;
int i,j,k;
int start, anzahl;
int index;
char c;
anzahl=derivedcomplete->firstfree;
start=0;
/* Die neuen Elemente der Liste befinden sich immer zwischen start und
   anzahl.
 */
while (start!=anzahl) 
   {
   for (i=start;i<anzahl;i++)
       {
       for (k=0;k<derivedgen->firstfree;k++)
           {
           hilf=(derived_TYP *)malloc(sizeof(derived_TYP));
	   hilf->nproduct = 0;
           hilf->element=mat_mul(derivedcomplete->list[i]->element,derivedgen->list[k]->element); 
          if (!is_in_list(derivedcomplete,hilf,&index))
             {
             hilf->nproduct=derivedcomplete->list[i]->nproduct+1;
             hilf->product=(int*)malloc(sizeof(int)*
                                hilf->nproduct);
             for (j=0;j<derivedcomplete->list[i]->nproduct;j++)
                 hilf->product[j]=derivedcomplete->list[i]->product[j];
             hilf->product[hilf->nproduct-1]=k;
             insert_list(derivedcomplete,hilf);   
             }
          else
             {
             free_derived(hilf);
             }
          }
       } 
   start=anzahl;
   anzahl=derivedcomplete->firstfree;
   }
}


/*--------------------------------------------------------------------------
   Funktion, das den Kommutator zweier Elemente berechnet (Spaltenkonvention)
   aba^-1b^-1.
   Benutzte Funktionen: (Matrixfunktionen)
      mat_mul
      mat_inv
      free_mat
      malloc
---------------------------------------------------------------------------*/
derived_TYP* kommutator(element1,element2)
matrix_TYP* element1, *element2;
{
derived_TYP* hilf;
matrix_TYP *hilf3, *hilf4, *hilf5;
hilf=(derived_TYP*) malloc(sizeof(derived_TYP));
hilf->nproduct=0;
hilf->product = NULL;
hilf->left = NIL;
hilf->right = NIL;
hilf3=mat_mul(element2, element1);
hilf4=mat_inv(hilf3);
hilf5=mat_mul(element1, element2);
hilf->element=mat_mul(hilf5,hilf4);
free_mat(hilf4);
free_mat(hilf3);
return (hilf);
}


/*-----------------------------------------------------------------------
   Die Prozedur derivedsg soll die Kommutatorgruppe einer durch Erzeuger
   definierten Gruppe bestimmen. Die Erzeuger der Untergruppe stehen in
   kommseriesgen[depth+1], die komplette Gruppe steht in 
   kommseriescomp[depth+1], wobei zu jedem Element die Zusammensetzung 
   als Produkt der Erzeuger mit abgespeichert ist. 
--------------------------------------------------------------------------*/
int derivedsg(kommseriesgen,depth, dimension,kommseriescomp)
derivedsg_TYP** kommseriesgen; 
int depth;
int dimension; 
derivedsg_TYP** kommseriescomp;
{
derived_TYP* kommnew;
int index;
int i, j;
int start;
int ende;
/* Initialisieren der Komplettliste (enthaelt die Einheitsmatrix)*/
kommseriescomp[depth+1]=init_derivedcomplete(dimension);
/* Initialisieren der Generatorliste (ist leer) */
kommseriesgen[depth+1]=init_derivedgen(dimension);
/* Erzeugnis aller aus den Erzeugern (der darueberliegenden Gruppe) 
   gebildeten Kommutatoren bilden. */
for (i=0;i<kommseriesgen[depth]->firstfree;i++)
    for (j=i+1;j<kommseriesgen[depth]->firstfree;j++)
    {
    kommnew=kommutator(kommseriesgen[depth]->list[i]->element,kommseriesgen[depth]->list[j]->element); 
    if (!is_in_list(kommseriescomp[depth+1],kommnew,&index)) 
       {
       insert_list(kommseriesgen[depth+1], kommnew);
       generate_sg(kommseriesgen[depth+1], kommseriescomp[depth+1]);
       }
    else
       {
       free_derived(kommnew);
       }
    }
start=0;
ende=kommseriesgen[depth+1]->firstfree;
/* Alle aus der vorherigen Schleife gewonnenen Generatoren werden mit
   den Generatoren der darueberliegenden Gruppe konjugiert. Das Ergebnis
   wird ggf. zum Erzeugendensystem hinzugenommen. Dieser Vorgang (i.e.
   konjugieren der Erzeuger) wird solange wiederholt, bis keine neuen
   Erzeuger entstehen. */
do
   {
   for (i=start;i<ende;i++)
      {
      for (j=0;j<kommseriesgen[depth]->firstfree;j++)
         {
         kommnew=konjugate(kommseriesgen[depth+1]->list[i],kommseriesgen[depth]->list[j]->element);
         if (!is_in_list(kommseriescomp[depth+1],kommnew,&index) )
            {
            insert_list(kommseriesgen[depth+1], kommnew);
            generate_sg(kommseriesgen[depth+1], kommseriescomp[depth+1]);
            }
         else
           free_derived(kommnew);
         }
      }
   start=ende;
   ende=kommseriesgen[depth+1]->firstfree;
}
while (start!=ende);
}
/*--------------------------------------------------------------------------
   Funktion zur Ausgabe der Relationen 
---------------------------------------------------------------------------*/
void dumprelators(relators, relator_no)
relator_TYP* relators;
int relator_no;
{
int i,j;
for (i=0;i<relator_no;i++)
    {
    printf("\n Relation no %d :\n  ",i);
    for (j=0;j<relators[i].lhsnproduct;j++)
        printf(" %d  ",relators[i].lhsproduct[j]);
    printf(" = ");
    for (j=0;j<relators[i].rhsnproduct;j++)
        printf(" %d  ",relators[i].rhsproduct[j]);
    printf("\n");
    }
printf("\n");
}

/* -------------------------------------------------------------------------
   Prozedur, die fuer einen neuen Generator (in der Praesentation), die
   Ordnung ueber der Darunterliegenden Gruppe bestimmt (genannt Relation
   vom Typ 1) und dann (Relationen vom Typ 2) ai*aj bestimmt fuer 
   Erzeuger ai,aj, i>j als Produkt
   von Erzeugern ak mit k<i und ai (ai ganz recht im Produkt).
   Dies ermoeglicht es alle Generatoren mit grossem Index nach rechts
   durchzuschieben.
--------------------------------------------------------------------------*/
void order(presentation, depth, genindex,  
           kommseriescomp)
presentation_TYP * presentation;
int depth;
int genindex; 
derivedsg_TYP** kommseriescomp;
{
derived_TYP* hilf;
matrix_TYP *hilf2;
int order_no;
int* product; 
int nproduct;
int index;
int i, j;
int k;
hilf=(derived_TYP*) malloc (sizeof(derived_TYP));
hilf->element=copy_mat(presentation->generators->list[presentation->
generators->firstfree-1]->element);
hilf->nproduct = 0;
hilf->product = NULL;
hilf->left = NIL;
hilf->right = NIL;
order_no=1;
do  {
    order_no++;
    hilf2=copy_mat(hilf->element);
    free_mat(hilf->element);
    hilf->element=mat_mul(hilf2,presentation->generators->list[presentation->generators->firstfree-1]->element);
    free_mat(hilf2);
    }
    while (!is_in_list(kommseriescomp[depth],hilf,&index));
free_derived(hilf);
if (presentation->norelators>=presentation->ext_factor*EXT_SIZE)
   {
   presentation->ext_factor++;
   presentation->relators = (relator_TYP*) realloc (presentation->relators, sizeof(relator_TYP)*EXT_SIZE*presentation->ext_factor);
   }
presentation->relators[presentation->norelators].lhsproduct = (int*) malloc (sizeof(int)*order_no);
for (i=0;i<order_no;i++)
    {
    presentation->relators[presentation->norelators].lhsproduct[i]=genindex;  
    }
product=(int*)malloc(sizeof(int)*kommseriescomp[depth]->list[index]->nproduct);
nproduct=kommseriescomp[depth]->list[index]->nproduct;
for (j=0;j<kommseriescomp[depth]->list[index]->nproduct;j++)
       product[j]=kommseriescomp[depth]->list[index]->product[j];
presentation->relators[presentation->norelators].lhsnproduct=order_no;
presentation->relators[presentation->norelators].rhsproduct=product;
presentation->relators[presentation->norelators].rhsnproduct=nproduct;
presentation->norelators++;
/* Bestimmen der Relationen "2. Art" */
k=presentation->norelators;
for (i=0;i<presentation->generators->firstfree-1;i++)
    {
    presentation->relators[k].lhsnproduct=2;
    presentation->relators[k].lhsproduct=(int*)malloc(sizeof(int)*2);
    presentation->relators[k].lhsproduct[0]=presentation->generators->firstfree-1;
    presentation->relators[k].lhsproduct[1]=i;
    hilf= konjugate(presentation->generators->list[i],presentation->generators->list[presentation->generators->firstfree-1]->element);
    is_in_list(kommseriescomp[depth],hilf,&index);
    presentation->relators[k].rhsnproduct=kommseriescomp[depth]->list[index]->nproduct+1;
    presentation->relators[k].rhsproduct=(int*)malloc(sizeof(int)*(presentation->relators[k].rhsnproduct+1));
    presentation->relators[k].rhsproduct[presentation->relators[k].rhsnproduct-1]=presentation->generators->firstfree-1;
    for (j=0;j<kommseriescomp[depth]->list[index]->nproduct;j++)
        presentation->relators[k].rhsproduct[j]=kommseriescomp[depth]->list[index]->product[j];
    k++;
    if (k>=presentation->ext_factor*EXT_SIZE)
       {
       presentation->ext_factor++;
       presentation->relators = (relator_TYP*) realloc (presentation->relators, sizeof(relator_TYP)*EXT_SIZE*presentation->ext_factor);
       }
    }
presentation->norelators=k;
}


/*--------------------------------------------------------------------------
   Funktion zur Ausgabe der Relationen auf der Datei presentation.dat
---------------------------------------------------------------------------*/
void dumprelators_file(outdat, relators, relator_no)
relator_TYP* relators;
int relator_no;
FILE* outdat;
{
int i,j;
fprintf(outdat,"\n%d\n",relator_no);
for (i=0;i<relator_no;i++)
    {
    for (j=0;j<relators[i].lhsnproduct;j++)
        fprintf(outdat," %d  ",relators[i].lhsproduct[j]+1);
    fprintf(outdat,"0\n");
    }
fprintf(outdat,"\n");
}

relator_TYP* init_relator()
{
relator_TYP* relations1;
relations1=(relator_TYP*)malloc(sizeof(relator_TYP)*EXT_SIZE);
return(relations1);
}

void dumpproduct(product, nproduct)
int * product;
int nproduct;
{
int i;
printf("\n");
for (i=0;i<nproduct;i++)
    printf(" %d  ",product[i]);
printf("\n");
}


/*-------------------------------------------------------------------------
  Prozedur, die fuer die in group abgelegten Gruppe, die nur durch Erzeuger
  gegeben ist, eine Praesentation bestimmt.
  Die Generatoren werden in generators zurueckgegeben.
  Die oben definierten Relationen 1.Art in relations1 und 2. Art in
  relations2. r2firstfree gibt die Anzahl der Relatoren 2. Art an.
  Groupgenold enthaelt die alten Erzeuger der Gruppe und ihre Darstellung
  als Produkt der neuen Erzeuger.
  Groupgennew enhaelt die neuen Erzeuger dargestellt als Produkt der alten.
---------------------------------------------------------------------------*/ 
int relations(group, presentation,
               groupgenold, groupgennew)
bravais_TYP* group;
presentation_TYP* presentation;
derivedsg_TYP *groupgenold;
derivedsg_TYP *groupgennew;
{
matrix_TYP* mattest;
derivedsg_TYP** kommseriesgen;
derivedsg_TYP **kommseriescomp;
derivedsg_TYP *groupcomplete;
derived_TYP* hilf;
int i,j,k;
int index;
char c;
int serieslength, generator_no;
presentation->norelators = 0;
/* Kommutatorreihe Bestimmen 
*/
kommseriesgen=(derivedsg_TYP **) malloc (sizeof(derivedsg_TYP*)*MAX_SERLEN);
kommseriescomp=(derivedsg_TYP **) malloc (sizeof(derivedsg_TYP*)*MAX_SERLEN);
kommseriesgen[0]=init_derivedgen(group->dim);
for (i=0;i<group->gen_no;i++)
    {
    hilf=(derived_TYP*) malloc(sizeof(derived_TYP));
    hilf->nproduct = 0;
    hilf->product = NULL;
    hilf->element=copy_mat(group->gen[i]);
    insert_list (kommseriesgen[0],hilf);
    hilf=(derived_TYP*) malloc(sizeof(derived_TYP));
    hilf->element=copy_mat(group->gen[i]);
    hilf->nproduct = 0;
    hilf->product = NULL;
    insert_list (groupgenold, hilf);
    }
j=0;
do   
  {
  derivedsg(kommseriesgen,j,group->dim, kommseriescomp); 
  fprintf(stderr," Derived subgroup %d , No of elements: %d, No of generators: %d\n",
      j+1, kommseriescomp[j+1]->firstfree, kommseriesgen[j+1]->firstfree);
  if ((j>0)&&(kommseriescomp[j+1]->firstfree==kommseriescomp[j]->firstfree))
     {
     printf(" Error: This group is not soluable!\n");
     return(-1);
     }
  j++;
  } while (kommseriesgen[j]->firstfree!=0);
serieslength=j;
generator_no=0;
/* Ende der Kommutatorreihenbestimmung
*/
/* Bestimmung der neuen Gruppengeneratoren und der Relationen, i.e. der
   Praesentation
*/
/* Ueber die gesamte Kommutatorreihe
*/
/* Speicherplatz wieder freigeben */
for (i=1;i<serieslength;i++)
    {
    free_derivedsgcomp(kommseriescomp[i]);
    }
for (i=serieslength;i>0;i--)
    {
    kommseriescomp[i-1]=copy_derivedsg(kommseriescomp[i]);
    free_derivedsg(kommseriesgen[i]);
    if (i<serieslength)
       free_derivedsgcomp(kommseriescomp[i+1]);
/* Ueber alle Generatoren dieses Kommutators
*/
    for (j=0;j<kommseriesgen[i-1]->firstfree;j++)
        {
        if (!is_in_list(kommseriescomp[i-1],kommseriesgen[i-1]->list[j],&index))
           {
           generator_no++;
           hilf=(derived_TYP*)malloc(sizeof(derived_TYP));
           hilf->element=copy_mat(kommseriesgen[i-1]->list[j]->element);
           hilf->nproduct = 0;
           hilf->product = NULL;
           insert_list(presentation->generators,hilf);
           hilf=(derived_TYP*)malloc(sizeof(derived_TYP));
           hilf->element=copy_mat(kommseriesgen[i-1]->list[j]->element);
           hilf->nproduct = 0;
           hilf->product = NULL;
           insert_list(groupgennew, hilf);
           order(presentation, 
                 i-1,
                 generator_no-1,
                 kommseriescomp); 
           generate_sg(presentation->generators,kommseriescomp[i-1]); 
           }
        }
    fprintf(stderr," no of generators %d , depth %d  \n",generator_no,i);
    }
/* Darstellung der urspruenglichen Erzeuger als Produkt der neuen Berechnen */
free_derivedsg(kommseriesgen[0]);
for (i=0;i<groupgenold->firstfree;i++)
    {
    is_in_list(kommseriescomp[0],groupgenold->list[i],&index);
    groupgenold->list[i]->product=copy_product(kommseriescomp[0]->list[index]);
    groupgenold->list[i]->nproduct=kommseriescomp[0]->list[index]->nproduct;
    }
/* Speicherplatz wieder freigeben */
for (i=0;i<2;i++)
    {
    free_derivedsgcomp(kommseriescomp[i]);
    }
free(kommseriesgen);
free(kommseriescomp);
groupcomplete=init_derivedcomplete(group->dim);
generate_sg(groupgenold, groupcomplete);
fprintf(stderr," Order of the Goup: %d\n",groupcomplete->firstfree);
/* Darstellung der neuen Erzeuger als Produkt der alten Berechnen
   d.h. Erzeugen der Gruppe mit den alten Erzeugern */
for (i=0;i<groupgennew->firstfree;i++)
    {
    is_in_list(groupcomplete, groupgennew->list[i],&index);
    groupgennew->list[i]->product=copy_product(groupcomplete->list[index]);
    groupgennew->list[i]->nproduct=groupcomplete->list[index]->nproduct;
    }
free_derivedsgcomp(groupcomplete);
return(0);
}



void shorter_product(ele_list)
derivedsg_TYP* ele_list;
{
int i,j;
for (i=0;i<ele_list->firstfree;i++)
    {
    if (ele_list->list[i]->product[0]==-1) 
       {
       for (j=0;j<ele_list->list[i]->nproduct-1;j++)
           ele_list->list[i]->product[j]=ele_list->list[i]->product[j+1];
       ele_list->list[i]->nproduct--;
       }
    }
}


main(argc, argv)
int argc;
char **argv;
{
bravais_TYP* testgroup;
int erg;
int i;
int found;
int dimension;
int* orders;
FILE* outdat;
presentation_TYP* presentation;
derivedsg_TYP* groupgenold, *groupgennew;
char file_out[150];

fprintf(stderr,"\n\n Program calculating a presentation of soluable groups \n");
fprintf(stderr,"\n (c) 1994 Lehrstuhl B Mathematik, RWTH Aachen\n\n\n");

/* inserted tilman 26/3/97 to conform with the rest of carat */
read_header(argc,argv);

if ((is_option('h') && optionnumber('h')) || FILEANZ < 1){
     printf("Usage\n ");
     printf(" Presentation file [-C] [-M] [-G]\n ");
     printf("\n ");
     printf("where file contains a bravais_TYP describing\n ");
     printf("a FINITE SOLUABLE group.\n ");
     printf("\n ");
     printf("Calculates a presentation of the group in file, and writes\n ");
     printf("it to the file <file>.carat.short .\n ");
     printf("\n ");
     printf("The options are:\n");
     printf("\n");
     printf(" -C   : Output the presentation in Cayley format on a file\n");
     printf("        called <file>.cayley.short .\n");
     printf(" -G   : Output the presentation in Gap format on a file\n");
     printf("        called <file>.gap.short .\n");
     printf(" -M   : Output the presentation in Magma format on a file\n");
     printf("        called <file>.magma.short .\n");
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

fprintf(stderr," The name of the input file: %s \n", FILENAMES[0]);

testgroup=get_bravais(FILENAMES[0]);
dimension = testgroup->dim;
presentation = (presentation_TYP*) malloc (sizeof(presentation));
presentation->relators=init_relator();
presentation->generators = init_derivedgen(testgroup->dim);
presentation->ext_factor = 1;
groupgenold=init_derivedgen(testgroup->dim);
groupgennew=init_derivedgen(testgroup->dim);
erg=relations(testgroup,presentation,
              groupgenold,groupgennew);
/*free_bravais(testgroup);*/
if (erg==0) 
   {
/*
   shorter_product(groupgenold);
   shorter_product(groupgennew);
   dumprelators(presentation->relators,presentation->norelators);
*/
   tietzetrans(groupgenold, groupgennew, presentation); 
/*
   printf(" Generators of the group \n");
   dumplist(presentation->generators);
*/
   /* this program isn't supported anymore
   sprintf(file_out,"%s.drin1",FILENAMES[0]);
   outdat=fopen(file_out,"a");
   dumplist_file(outdat, presentation->generators);
   fclose(outdat); */
/*
   dumprelators(presentation->relators,presentation->norelators);
*/
   /* this program isn't supported anymore 
   sprintf(file_out,"%s.rel",FILENAMES[0]);
   outdat=fopen(file_out,"a");
   dumprelators_file(outdat, presentation->relators,presentation->norelators);
   fclose(outdat); */

   if (is_option('C')){
      sprintf(file_out,"%s.cayley",FILENAMES[0]);
      cayley_out(presentation,file_out);
   }

   orders = (int*) malloc (sizeof(int)*presentation->generators->firstfree);
   calc_orders(presentation, orders, dimension);
   shorten_presentation(presentation, orders);

   /* this program is not supported anymore
   sprintf(file_out,"%s.rel.short",FILENAMES[0]);
   outdat=fopen(file_out,"a");
   dumprelators_file(outdat, presentation->relators,presentation->norelators);
   fclose(outdat); */

   /* reorganized the output tilman 26/3/97 */
   /* sprintf(file_out,"%s.carat.short",FILENAMES[0]);
   fprintf(stderr,"\n Output written on: %s \n", file_out); */
   carat_out(presentation,NULL);

   if (is_option('M')){
      sprintf(file_out,"%s.magma.short",FILENAMES[0]);
      fprintf(stderr,"\n Output written on: %s \n", file_out);
      magma_out(presentation,file_out);
   }

   if (is_option('G')){
      sprintf(file_out,"%s.gap.short",FILENAMES[0]);
      fprintf(stderr,"\n Output written on: %s \n", file_out);
      gap_out(presentation,file_out);
   }

   if (is_option('C')){
      sprintf(file_out,"%s.cayley.short",FILENAMES[0]);
      fprintf(stderr,"\n Output written on: %s \n", file_out);
      cayley_out(presentation,file_out);
   }

   }
else{
   /* the group is not soluble, and it has been told to the audience
      already */
  exit(3);
}

free_derivedsgcomp(groupgenold);
free_derivedsgcomp(groupgennew);
free_derivedsg(presentation->generators);
for (i=0;i<presentation->norelators;i++)
    {
    if (presentation->relators[i].lhsproduct!=NULL)
       free(presentation->relators[i].lhsproduct);
    if (presentation->relators[i].rhsproduct!=NULL &&
        presentation->relators[i].rhsnproduct != 0)
       free(presentation->relators[i].rhsproduct);
    }
free(presentation->relators);


exit(0);
}
