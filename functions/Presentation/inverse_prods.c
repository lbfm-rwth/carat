#include "typedef.h"
/* #include "const.h" */
#include "dertype.h"
#include "tietzetrans.h"

/*--------------------------------------------------------------------------
   Prozedur, die zu jedem Generator in presentation das Inverse berechnet,
   als Wort ebendieses Generators. Also a^(-1) = a^(n-1), wobei a der 
   Generator ist und n die Ordnung deselben. Daraufhin wird mit diesen
   Daten die Prozedur shorter (siehe dazu shorter.c) aufgerufen.
   Parameter:
      presentation - Praesentation
      orders - Ordnung der Generatoren
--------------------------------------------------------------------------*/ 
void inverse_prods(presentation, orders)
presentation_TYP* presentation;
int* orders;
{
int** inverse_product;
int * inverse_nproduct;
int i;
int j;
inverse_product = (int**) malloc (sizeof(int*)*presentation->generators->firstfree);
inverse_nproduct = (int*) malloc (sizeof(int)*presentation->generators->firstfree);
for (i=0;i<presentation->generators->firstfree;i++)
    {
    inverse_product[i]=(int*) malloc (sizeof(int)*(orders[i]-1));
    inverse_nproduct[i]=orders[i]-1;
    for (j=0;j<orders[i]-1;j++)
        {
        inverse_product[i][j]=i;
        }
    }
shorter(presentation, inverse_product, inverse_nproduct);
for (i=0;i<presentation->generators->firstfree;i++)
    {
    free(inverse_product[i]);
    }
free(inverse_product);
free(inverse_nproduct);
}
