
#include "typedef.h"

#include "dertype.h"
#include "tietzetrans.h"
#include "roundtype.h"
#define maxgen 9

/*-------------------------------------------------------------------------
   Prozedur zur Ausgabe einer Praesentation im GAP-Format.
   Parameter:
        presentation - Praesentation
        filename - Dateiname fuer die Ausgabe
---------------------------------------------------------------------------*/
void gap_out (presentation,filename)
presentation_TYP* presentation;
char* filename;
{
char* gens;
int i,j;
FILE* outfile;
outfile = fopen(filename,"w");
gens=(char*) malloc (sizeof(char)*maxgen);
gens[0]='1';
gens[1]='2';
gens[2]='3';
gens[3]='4';
gens[4]='5';
gens[5]='6';
gens[6]='7';
gens[7]='8';
gens[8]='9';
fprintf(outfile," G := FreeGroup(%d);\n", presentation->generators->firstfree);
fprintf(outfile," H := G/[");
for (i=0;i<presentation->norelators;i++)
    {
    if (presentation->relators[i].lhsnproduct>0)
       fprintf(outfile,"G.%c",gens[presentation->relators[i].lhsproduct[0]]);
    for (j=1;j<presentation->relators[i].lhsnproduct;j++)
        {
        fprintf(outfile,"*G.%c",gens[presentation->relators[i].lhsproduct[j]]);
        if ((j*2)%70 == 0) fprintf(outfile, "\n");
        }
    if (presentation->relators[i].rhsnproduct>0)
       {
       fprintf(outfile," = ");
       fprintf(outfile,"G.%c",gens[presentation->relators[i].rhsproduct[0]]);
       }
    for (j=1;j<presentation->relators[i].rhsnproduct;j++)
        {
        fprintf(outfile,"*G.%c",gens[presentation->relators[i].rhsproduct[j]]);
        if ((j*2)%70 == 0) fprintf(outfile, "\n");
        }
    if (i<presentation->norelators-1)
       {
       fprintf(outfile,",\n");
       }
    }
fprintf(outfile,"];\n");
fprintf(outfile,"Size(H);");
fclose(outfile);
free(gens);
} 

