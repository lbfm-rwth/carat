#include "typedef.h"
/* #include "const.h" */
#include "dertype.h"
#include "tietzetrans.h"

/*-------------------------------------------------------------------------
   Prozedur, die in den Relatoren ueberprueft, ob Bloecke gleicher
   Generatoren laenger sind, als die Ordnung dieses Elements und ggf.
   die ueberfluessigen Elemente des Relators loescht.
   (dies war mal ein Teil von Tietzetrans.c, uwrde aber aufgrund von
    Mehrfachverwendung extrahiert).    
    Parameter: 
       presentation - Praesentation
       inverse_product - Liste der Darstellungen der Inversen der 
                 Generatoren.
       inverse_nproduct - Laenge der Darstellungen in inverse_product.
--------------------------------------------------------------------------*/
void shorter(presentation, inverse_product, inverse_nproduct)
presentation_TYP* presentation;
int* inverse_nproduct;
int** inverse_product;
{
int i;
int blockcount;
int lauf;
int number;
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
 /* in den Ordnungsrelationen muss natuerlich der Block einmal stehen
    bleiben. */
          if (blockcount==presentation->relators[i].lhsnproduct)
             {
             number-=(inverse_nproduct[presentation->relators[i].lhsproduct[lauf]]+1);
             }
          number-=(inverse_nproduct[presentation->relators[i].lhsproduct[lauf]]+1);
          remove_product(&presentation->relators[i].lhsproduct,
                &presentation->relators[i].lhsnproduct,lauf,number);
          lauf=0;
          }
       else
          lauf+=blockcount;
       }
    }
}
