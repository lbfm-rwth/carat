#include "typedef.h"

/*-----------------------------------------------------------------------
   Funktion, die angibt, ob ein File existiert oder nicht.
   Parameter:
        filename - Name der Datei
   Rueckgabewert: 
        1 - Datei existiert
        0 - Datei existiert nicht
------------------------------------------------------------------------*/
int exist(filename)
char *filename;
{
FILE* tf;
if (( tf=fopen(filename, "r"))==NULL)
   {
   return(0);
   }
fclose(tf);
return(1);
}
