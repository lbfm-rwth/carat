
#include "typedef.h"

#include "dertype.h"
#include "tietzetrans.h"
#include "roundtype.h"
#define maxgen 26

/*--------------------------------------------------------------------------
  Prozedur zur Ausgabe einer Praesentation im Cayley-Format.
  Parameter:
      presentation - Praesentation
      filename - Name des outputfiles
---------------------------------------------------------------------------*/
void cayley_out (presentation,filename)
presentation_TYP* presentation;
char* filename;
{
char* gens;
int i,j;
FILE* outfile;
outfile = fopen(filename,"w");
gens=(char*) malloc (sizeof(char)*maxgen);
gens[0]='a';
gens[1]='b';
gens[2]='z';
gens[3]='c';
gens[4]='d';
gens[5]='e';
gens[6]='f';
gens[7]='y';
gens[8]='h';
gens[9]='i';
gens[10]='j';
gens[11]='k';
gens[12]='l';
gens[13]='m';
gens[14]='n';
gens[15]='o';
gens[16]='p';
gens[17]='q';
gens[18]='r';
gens[19]='s';
gens[20]='y';
gens[21]='t';
gens[22]='u';
gens[23]='v';
gens[24]='w';
gens[25]='x';
fprintf(outfile," G = free(a");
for (i=1;i<presentation->generators->firstfree;i++)
    {
    fprintf(outfile,",%c",gens[i]);
    }
fprintf(outfile,");\n G.relations: \n");
for (i=0;i<presentation->norelators;i++)
    {
    if (presentation->relators[i].lhsnproduct>0) 
      fprintf(outfile,"%c",gens[presentation->relators[i].lhsproduct[0]]);
    for (j=1;j<presentation->relators[i].lhsnproduct;j++)
        {
        fprintf(outfile,"*%c",gens[presentation->relators[i].lhsproduct[j]]);
        if ((j*2)%70 == 0) fprintf(outfile, "\n");
        }
    if (presentation->relators[i].rhsnproduct>0)
       {
       fprintf(outfile," = ");
       fprintf(outfile,"%c",gens[presentation->relators[i].rhsproduct[0]]);
       }
    for (j=1;j<presentation->relators[i].rhsnproduct;j++)
        {
        fprintf(outfile,"*%c",gens[presentation->relators[i].rhsproduct[j]]);
        if ((j*2)%70 == 0) fprintf(outfile, "\n");
        }
    if (i<presentation->norelators-1)
       {
       fprintf(outfile,",\n");
       }
    else
       {
       fprintf(outfile,";\n");
       }
    }
fprintf(outfile," print order(G);\n");
fclose(outfile);
free(gens);
} 

