
#include "typedef.h"

#include "dertype.h"
#include "tietzetrans.h"
#include "roundtype.h"
#define maxgen 9

/*-------------------------------------------------------------------------
   Prozedur zur Ausgabe einer Praesentation im CARAT-Format.
   Parameter:
        presentation - Praesentation
        filename - Dateiname fuer die Ausgabe
---------------------------------------------------------------------------*/
void carat_out (presentation,filename)
presentation_TYP* presentation;
char* filename;
{

extern matrix_TYP* init_mat(int, int,char*);
char* gens;
int i,j;
int max_cols;
matrix_TYP* erg;
FILE* outfile;

/* changed to stand NULL for stdout */
if (filename == NULL)
   outfile = stdout;
else
   outfile = fopen(filename,"w");

max_cols = 0;
for (i=0;i<presentation->norelators;i++)
    {
    if (presentation->relators[i].lhsnproduct>max_cols)
	max_cols = presentation->relators[i].lhsnproduct;
    }
erg = init_mat(presentation->norelators, max_cols,"k");
for (i=0;i<presentation->norelators;i++)
    	{
	for (j=0;j<max_cols;j++)
		{
		erg->array.SZ[i][j] = 0;
		}
	}
for (i=0;i<presentation->norelators;i++)
    	{
	for (j=0;j<presentation->relators[i].lhsnproduct;j++)
		{
		erg->array.SZ[i][j] = 1+presentation->relators[i].lhsproduct[j];
		}
	}	
dumpmat_file(outfile, erg);

/* changed to deal with stdout */
if (filename != NULL) fclose(outfile);

free_mat(erg);
} 

