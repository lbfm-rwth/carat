#include "typedef.h"
#include "dertype.h"

/*-------------------------------------------------------------------------
   Funktion, die eine Ordnung auf Matrizen realisiert. Die beiden
   Elemente werden verglichen und wenn element1 groesser (in Bezug auf
   diese Ordnung) ist als element2, wird 1 zurueckgegeben, sonst 0. 
-------------------------------------------------------------------------*/
int is_equal (element1,element2)
derived_TYP* element1,* element2;
{
int i,j;
for (i=0;i<element1->element->rows;i++)
    for (j=0;j<element1->element->cols;j++)
        if (element2->element->array.SZ[i][j] !=
              element1->element->array.SZ[i][j])
           {
           return (FALSE);
           }
return(TRUE);
}
