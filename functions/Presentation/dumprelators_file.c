#include "typedef.h"
/* #include "const.h" */
#include "dertype.h"
#include "tietzetrans.h"
#define MAX_SERLEN 10


/*--------------------------------------------------------------------------
   Funktion zur Ausgabe der Relationen auf der Datei outdat 
   Parameter:
        outdat - Pointer auf den Ausgabefile
        relators - Relatoren
        relator_no - Anzahl Relatore
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
