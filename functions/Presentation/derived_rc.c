#include "typedef.h"
#include "dertype.h"
#include <stdio.h>

extern matrix_TYP *init_mat();
extern int mat_vergleich();
extern matrix_TYP *copy_mat();

/* Prototypenliste -------------------------------------------------------*/
derivedsg_TYP* init_derivedcomplete (int);
derivedsg_TYP* init_derivedgen (int);
derived_TYP* copy_derived(derived_TYP*);
int* copy_product(derived_TYP*);
derivedsg_TYP* copy_derivedsg(derivedsg_TYP*);
void dumpmat_file(FILE*, matrix_TYP*);
void dumpmat(matrix_TYP*);
void dumplist_file(FILE* , derivedsg_TYP*);
void dumplist(derivedsg_TYP*);
void insert_list (derivedsg_TYP*, derived_TYP*);
void free_derivedsg(derivedsg_TYP*);
void free_derivedsgcomp(derivedsg_TYP*);
void free_derived(derived_TYP*);
/*
int isequal(derived_TYP*, derived_TYP*);
*/
int isgreater (derived_TYP*,derived_TYP*);
int is_in_list (derivedsg_TYP*, derived_TYP*, int*);
void delete_ele (derivedsg_TYP*, derived_TYP*, int*);

/*------------------------------------------------------------------------
   Initialisieren einer Liste, in der alle Elemente einer Untergruppe
   abgespeichert werden. Dabei wird das 1-Element der Gruppe bereits
   eingefuegt. 
-------------------------------------------------------------------------*/
derivedsg_TYP* init_derivedcomplete (dimension)
int dimension;
{
derivedsg_TYP* hilf;
int i;
hilf=(derivedsg_TYP *) malloc (sizeof(derivedsg_TYP));
hilf->list=(derived_TYP **) malloc (sizeof(derived_TYP*)*EXT_SIZE);
hilf->sizemult=1;
hilf->list[0]=(derived_TYP*)malloc(sizeof(derived_TYP));
hilf->list[0]->element=init_mat(dimension,dimension,"k");
hilf->list[0]->right=NIL;
hilf->list[0]->left=NIL;
hilf->firstfree=1;
for (i=0;i<dimension;i++)
         hilf->list[0]->element->array.SZ[i][i]=1;
hilf->list[0]->product=(int*)malloc(sizeof(int));
hilf->list[0]->product[0]=NIL;
hilf->list[0]->nproduct=0;
return (hilf);
}

/* ------------------------------------------------------------------------
   Initialisieren einer Generatorliste. Hier wird nur Speicherplatz fuer
   die Zeiger reserviert. Die Liste wird leer zurueckgegeben. 
---------------------------------------------------------------------------*/
derivedsg_TYP* init_derivedgen (dimension)
int dimension;
{
derivedsg_TYP* hilf;
hilf=(derivedsg_TYP *) malloc (sizeof(derivedsg_TYP));
hilf->list=(derived_TYP **)malloc(sizeof(derived_TYP*)*EXT_SIZE);
hilf->firstfree=0;
hilf->sizemult=1;
return(hilf);
}

/*------------------------------------------------------------------------- 
   Funktion, die ein Gruppenelement kopiert. Der Inhalt des uebergebenen
   Elements wird dupliziert und ein Zeiger auf dea Duplikat zurueckgegeben.
---------------------------------------------------------------------------*/
derived_TYP* copy_derived(oriderived)
derived_TYP* oriderived;
{
int i;
derived_TYP* hilf;
hilf=(derived_TYP*)malloc(sizeof(derived_TYP));
hilf->element=copy_mat(oriderived->element);
hilf->nproduct=oriderived->nproduct;
hilf->product=(int*) malloc (hilf->nproduct*sizeof(int));
hilf->left=oriderived->left;
hilf->right=oriderived->right;
for (i=0;i<hilf->nproduct;i++)
    hilf->product[i]=oriderived->product[i];
return(hilf);
}

/*---------------------------------------------------------------------------
   Funktion, die die Darstellung eines Elements als Produkt der Erzeuger
   kopiert.
----------------------------------------------------------------------------*/
int* copy_product(oriderived)
derived_TYP* oriderived;
{
int* hilf;
int i;
hilf=(int*) malloc (oriderived->nproduct*sizeof(int));
for (i=0;i<oriderived->nproduct;i++)
    hilf[i]=oriderived->product[i];
return(hilf);
}

/*---------------------------------------------------------------------------
   Funktion um eine Liste zu kopieren. 
----------------------------------------------------------------------------*/
derivedsg_TYP* copy_derivedsg(oriderivedsg)
derivedsg_TYP* oriderivedsg;
{
int i;
derivedsg_TYP* copysg;
copysg=(derivedsg_TYP *) malloc (sizeof(derivedsg_TYP));
copysg->list=(derived_TYP **)malloc(sizeof(derived_TYP*)*EXT_SIZE*oriderivedsg->sizemult);
copysg->firstfree=oriderivedsg->firstfree;
copysg->sizemult=oriderivedsg->sizemult;
for (i=0;i<oriderivedsg->firstfree;i++)
    {
    copysg->list[i]=copy_derived(oriderivedsg->list[i]); 
    }
return(copysg);
}



/*--------------------------------------------------------------------------
   Hilfsprozeduren um eine Matrix auf Datei auszugeben.
----------------------------------------------------------------------------*/
void dumpmat_file(outdat, matrix)
FILE* outdat;
matrix_TYP* matrix;
{
int i,j;
fprintf(outdat, "%dx%d\n",matrix->rows, matrix->cols);
for (i=0;i<matrix->rows;i++)
 	{
	for (j=0;j<matrix->cols;j++)
		{
		fprintf(outdat, "%d ",matrix->array.SZ[i][j]);
		}
        fprintf(outdat, "\n");
	}
}
/*--------------------------------------------------------------------------
   Hilfsprozeduren um eine Matrix auszugeben.
----------------------------------------------------------------------------*/
void dumpmat(matrix)
matrix_TYP* matrix;
{
int i,j;
for (i=0;i<matrix->rows;i++)
 	{
	for (j=0;j<matrix->cols;j++)
		{
		printf("%d ",matrix->array.SZ[i][j]);
		}
        printf("\n");
	}
}

/* ------------------------------------------------------------------------
   Hilfsprozedur um eine Liste (d.h. die darin gespeicherten Matrizen)
   auszugeben (auf einer Datei). 
-------------------------------------------------------------------------*/
void dumplist_file(outdat, derivedsg)
FILE* outdat;
derivedsg_TYP *derivedsg;
{
int i;
char c;
fprintf(outdat,"\n%d\n", derivedsg->list[0]->element->rows);
fprintf(outdat,"%d\n", derivedsg->firstfree);
for (i=0;i<derivedsg->firstfree;i++)
	{
        fprintf(outdat," \n");
	dumpmat_file(outdat, derivedsg->list[i]->element);
	}
}
/* ------------------------------------------------------------------------
   Hilfsprozedur um eine Liste (d.h. die darin gespeicherten Matrizen)
   auszugeben. 
-------------------------------------------------------------------------*/
void dumplist(derivedsg)
derivedsg_TYP *derivedsg;
{
int i;
char c;
for (i=0;i<derivedsg->firstfree;i++)
	{
        printf(" \n");
	dumpmat(derivedsg->list[i]->element);
         printf(" left= %d, right= %d",derivedsg->list[i]->left,derivedsg->list[i]->right);
       /* c=getchar();*/
	}
}

/*------------------------------------------------------------------------ 
   Speicherplatz, der von einer Generatorliste belegt wurde wird wieder
   freigegeben. 
-------------------------------------------------------------------------*/
void free_derivedsg(derivedsg)
derivedsg_TYP* derivedsg;
{
int i;
for (i=0;i<derivedsg->firstfree;i++)
    {
    free_mat(derivedsg->list[i]->element);
    if (derivedsg->list[i]->nproduct>0)
       free(derivedsg->list[i]->product);
    free(derivedsg->list[i]);
    }
free(derivedsg->list);
free(derivedsg);
}

/*------------------------------------------------------------------------ 
   Speicherplatz, der von einer Elementliste belegt wurde wird wieder
   freigegeben.  
-------------------------------------------------------------------------*/
void free_derivedsgcomp(derivedsg)
derivedsg_TYP* derivedsg;
{
int i;
for (i=0;i<derivedsg->firstfree;i++)
    {
    free_mat(derivedsg->list[i]->element);
    free(derivedsg->list[i]->product);
    free(derivedsg->list[i]);
    }
free(derivedsg->list);
free(derivedsg);
}

/*-------------------------------------------------------------------------
   Speicherplatz eines Elements wieder freigeben 
-------------------------------------------------------------------------*/
void free_derived(derivedele)
derived_TYP* derivedele;
{
free_mat(derivedele->element);
if (derivedele->nproduct!=0)
	free(derivedele->product);
free(derivedele);
}

/*------------------------------------------------------------------------- 
  Vergleich zweier Elemente, dabei wird die Funktion mat_vergleich
  genutzt, die die Matrizen, der Elemente vergleicht. 
-------------------------------------------------------------------------*/
/******************************
int isequal(element1, element2)
derived_TYP* element1,* element2;
{
return(!mat_vergleich(element1->element, element2->element));
}
*********************************/

/*-------------------------------------------------------------------------
   Funktion, die eine Ordnung auf Matrizen realisiert. Die beiden
   Elemente werden verglichen und wenn element1 groesser (in Bezug auf
   diese Ordnung) ist als element2, wird 1 zurueckgegeben, sonst 0. 
-------------------------------------------------------------------------*/
int isgreater (element1,element2)
derived_TYP* element1,* element2;
{
int i,j;
for (i=0;i<element1->element->rows;i++)
    for (j=0;j<element1->element->cols;j++)
        if (element2->element->array.SZ[i][j] !=
              element1->element->array.SZ[i][j])
           {
            if (element2->element->array.SZ[i][j] <
                 element1->element->array.SZ[i][j])
              return (TRUE);
           else return (FALSE);
           }
}

/*-----------------------------------------------------------------------
   Funktion, die Ueberprueft, ob element in derivedlist enthalten ist.
   Falls ja, wird 1 zurueckgegeben und in index die Position in der 
   Liste auf der das Element zu finden ist abgelegt. 
   Linear suchende Version.
----------------------------------------------------------------------------*/
int is_in_list_l (derivedlist, element, index)
derivedsg_TYP* derivedlist;
derived_TYP* element;
int* index;
{
int i;
for (i=0;i<derivedlist->firstfree;i++)
   {
   if (is_equal(derivedlist->list[i],element))
      {
      *index=i;
      return(TRUE);
      }
   }
return(FALSE);
}


/*-----------------------------------------------------------------------
   Funktion, die Ueberprueft, ob element in derivedlist enthalten ist.
   Falls ja, wird 1 zurueckgegeben und in index die Position in der 
   Liste auf der das Element zu finden ist abgelegt.
----------------------------------------------------------------------------*/
int is_in_list (derivedlist, element, index)
derivedsg_TYP* derivedlist;
derived_TYP* element;
int* index;
{
int lauf=0;
if (derivedlist->firstfree==0) return (FALSE);
while (lauf!=NIL) 
   {
   if (is_equal(derivedlist->list[lauf],element))
      {
      *index=lauf;
      return(TRUE);
      }
   if (isgreater(derivedlist->list[lauf],element))
      lauf=derivedlist->list[lauf]->right;
   else
      lauf=derivedlist->list[lauf]->left;
   }
return(FALSE);
}


void delete_ele (derivedlist, element, index)
derivedsg_TYP* derivedlist;
derived_TYP* element;
int* index;
{
int i;
int lauf=0;
int laufalt=0, laufuralt=0;
matrix_TYP* hilf;
derived_TYP* hilfderived;
while (lauf!=NIL) 
   {
   if (is_equal(derivedlist->list[lauf],element))
      {
      hilf=derivedlist->list[lauf]->element;
      laufuralt=laufalt;
      laufalt=lauf;
      lauf=derivedlist->list[laufalt]->right;
      while (lauf!=NIL)
            {
            derivedlist->list[laufalt]->element=derivedlist->list[lauf]->element;
            laufuralt=laufalt;
            laufalt=lauf;
            lauf=derivedlist->list[laufalt]->right; 
            }
      if (laufuralt!=laufalt)
         {
         if (derivedlist->list[laufalt]->left!=NIL) 
            {
            derivedlist->list[laufuralt]->right=derivedlist->list[laufalt]->left;
            }
         else
            {
            derivedlist->list[laufuralt]->right=NIL;
            }
         derivedlist->list[laufalt]->element=hilf;
         hilfderived=derivedlist->list[laufalt];
         *index=laufalt;
         for (i=laufalt;i<derivedlist->firstfree-1;i++)
             derivedlist->list[i]=derivedlist->list[i+1];
         derivedlist->list[derivedlist->firstfree-1]=hilfderived;
         free_derived(derivedlist->list[derivedlist->firstfree-1]);
         }
      else
         {
/* Das zu loeschende Element hat keinen rechten Sohn und keinen Vater, also
   kann es ohne Ruecksicht auf die Baumstruktur aus dem Feld geloescht werden
*/
         hilfderived=derivedlist->list[laufalt];
         *index=laufalt;
         for (i=laufalt;i<derivedlist->firstfree-1;i++)
             derivedlist->list[i]=derivedlist->list[i+1];
         derivedlist->list[derivedlist->firstfree-1]=hilfderived;
         free_derived(derivedlist->list[derivedlist->firstfree-1]);
         }
      derivedlist->firstfree--;
      return;
      }
   laufuralt=laufalt;
   laufalt=lauf;
   if (isgreater(derivedlist->list[lauf],element))
      lauf=derivedlist->list[lauf]->right;
   else
      lauf=derivedlist->list[lauf]->left;
   }
return;
}

/*----------------------------------------------------------------------------
   Funktion, die element in derivedlist einfuegt. Eine Ueberpruefung, ob
   das Element bereits in der Liste ist wird nicht gemacht. Das Element
   wird an der Postition firstfree eingefuegt, und die Indexverpointerung
   des "Baumes" aktualisiert.
-----------------------------------------------------------------------------*/
void insert_list (derivedlist, element)
derivedsg_TYP* derivedlist;
derived_TYP* element;
{
int lauf, laufalt;
if (derivedlist->firstfree >= derivedlist->sizemult*EXT_SIZE)
  {
  derivedlist->sizemult++;
  derivedlist->list=(derived_TYP**)
                    realloc(derivedlist->list,derivedlist->sizemult*
                            sizeof(derived_TYP*)
                            *EXT_SIZE);
  }
derivedlist->list[derivedlist->firstfree]=element;
derivedlist->list[derivedlist->firstfree]->left=NIL;
derivedlist->list[derivedlist->firstfree]->right=NIL;
derivedlist->firstfree++;
if (derivedlist->firstfree!=1) 
   {
   lauf=0;
   laufalt=0;
   while (lauf!=NIL) 
      {
      laufalt=lauf;
      if (isgreater(derivedlist->list[lauf],element))
         lauf=derivedlist->list[lauf]->right;
      else
         lauf=derivedlist->list[lauf]->left;
      }
   if (isgreater(derivedlist->list[laufalt],element))
      derivedlist->list[laufalt]->right=derivedlist->firstfree-1;
   else
      derivedlist->list[laufalt]->left=derivedlist->firstfree-1;
   }
/*
printf("\n in insert_list \n");
dumplist(derivedlist);
*/
}
