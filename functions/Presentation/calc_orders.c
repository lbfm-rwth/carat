#include "typedef.h"
#include "dertype.h"
#include "tietzetrans.h"

extern int is_identity( presentation_TYP*, int);

/*-------------------------------------------------------------------------
   Programm zur Berechnung der Ordnungen der Erzeuger in presentation,
   also der Elemente in presentation->generators.
   Dabei wird fuer jeden Erzeuger der Gruppe ein Relator hinzugefuegt
   (in der Form a^n, wobei a der Erzeuger ist und n die Ordnung des 
   Erzeugers.
   Parameter:
       presentation - Praesentation (Ein/Ausgabe)
       orders - Ordnungen der Generatoren (Ausgabe)
       dim - Dimension des Raumes (Eingabe)
--------------------------------------------------------------------------*/
void calc_orders( presentation,orders, dim )
presentation_TYP* presentation;
int* orders;
int dim;
{
int i;
int k;
for (i=0; i<presentation->generators->firstfree;i++)
    {
    if (presentation->norelators>=(presentation->ext_factor*EXT_SIZE))
       {
       presentation->ext_factor++;
       presentation->relators=(relator_TYP*) realloc (presentation->relators,
          sizeof(relator_TYP)*presentation->ext_factor*EXT_SIZE);
       }
   /* Schreibabkuerzung */
    k=presentation->norelators;
    presentation->relators[k].lhsnproduct = 0;
    presentation->relators[k].lhsproduct = (int*) malloc (sizeof(int));
    presentation->relators[k].rhsnproduct = 0;
    do
       { 
       presentation->relators[k].lhsnproduct++;
       presentation->relators[k].lhsproduct = (int*) realloc 
            (presentation->relators[k].lhsproduct, sizeof(int)*
               presentation->relators[k].lhsnproduct);
       presentation->relators[k].lhsproduct[
           presentation->relators[k].lhsnproduct-1]
           = i;
       } 
    while (!is_identity( presentation, dim)); 
    orders[i] = presentation->relators[k].lhsnproduct;
    presentation->norelators++;
    }
}

/*-------------------------------------------------------------------------
    Prozedur zur Berechnung der Ordnung eines der Erzeuger in der 
    Praesentation und zwar des letzten in der Liste presentation->generators
    Parameter:
        presentation - Praesentation (Ein/Ausgabe)
        orders - Ordnungen der Erzeuger (Ausgabe)
        dim - Dimension des Raumes (Eingabe)
-------------------------------------------------------------------------*/
void calc_orders_s( presentation,orders, dim )
presentation_TYP* presentation;
int* orders;
int dim;
{
int i;
int k;
i=presentation->generators->firstfree-1;
if (presentation->norelators>=(presentation->ext_factor*EXT_SIZE))
   {
   presentation->ext_factor++;
   presentation->relators=(relator_TYP*) realloc (presentation->relators,
      sizeof(relator_TYP)*presentation->ext_factor*EXT_SIZE);
   }
/* Schreibabkuerzung */
k=presentation->norelators;
presentation->relators[k].lhsnproduct = 0;
presentation->relators[k].lhsproduct = (int*) malloc (sizeof(int));
presentation->relators[k].rhsnproduct = 0;
do
   { 
   presentation->relators[k].lhsnproduct++;
   presentation->relators[k].lhsproduct = (int*) realloc 
        (presentation->relators[k].lhsproduct, sizeof(int)*
           presentation->relators[k].lhsnproduct);
   presentation->relators[k].lhsproduct[
       presentation->relators[k].lhsnproduct-1]
       = i;
   } 
while (!is_identity( presentation, dim)); 
orders[i] = presentation->relators[k].lhsnproduct;
presentation->norelators++;
}
