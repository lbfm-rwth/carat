#include "typedef.h"
#include "dertype.h"
#include "tietzetrans.h"

/*--------------------------------------------------------------------------
   Funktion zur Ausgabe der Relationen 
   Parameter:
       relators - Relatoren
       relator_no - Anzahl Relatoren
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

