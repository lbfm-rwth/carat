#include "typedef.h"
/* #include "const.h" */
#include "dertype.h"
#include "tietzetrans.h"

/*---------------------------------------------------------------------------
   Diese Prozedur ersetzt in der gesamten Praesentation das auftauchen des
   Generators gen durch die rechte Seite des Relators relatorind. Dieser
   Relator sollte ein Darstellung des Generators als Produkt von anderen
   Generatoren beinhalten. Dies kann sichergestellt werden, durch vor-
   heriges aufrufen von which_relator.
---------------------------------------------------------------------------*/
void substitute_inverse (presentation, generator, order)
presentation_TYP* presentation;
int generator;
int order;
{
int i, j;
int * inverse_product;
int inverse_nproduct;
inverse_nproduct = order -1;
inverse_product = (int*) malloc (sizeof(int)*inverse_nproduct);
for (i=0;i<inverse_nproduct;i++)
    {
    inverse_product[i] = generator;
    }
for (i=0;i<presentation->norelators;i++)
    {
    j=0;
    while (j<presentation->relators[i].lhsnproduct)
        {
        if (presentation->relators[i].lhsproduct[j]==generator)
           {
           insert_in_product(&presentation->relators[i].lhsproduct, j, &presentation->relators[i].lhsnproduct,inverse_product, inverse_nproduct);
           j+=inverse_nproduct;
           }
        else
           {
            j++;
           }
        }
    j=0;
    while (j<presentation->relators[i].rhsnproduct)
        {
        if (presentation->relators[i].rhsproduct[j]==generator)
           {
           insert_in_product(&presentation->relators[i].rhsproduct, j, &presentation->relators[i].rhsnproduct, inverse_product, inverse_nproduct);
           j+= inverse_nproduct;
           }
        else
           {
           j++;
           }
        }
    }
}

