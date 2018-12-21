#include "typedef.h"
#include "tools.h"
#include"getput.h"
#include"datei.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  get_zentr.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ void get_zentr(B)
@ symbol_out *B;
@
@    Another function for the 'Datei' programm.
@    It reads matrices form the file 'B->fn' and 
@    stores them in 'B->zentr'.
@    Reads additionally a filename stated beyond the matrices
@    and stores it in B->fn.
@---------------------------------------------------------------------------
@
\**************************************************************************/
void get_zentr(B)
symbol_out *B;
{  

char *file_name;
boolean header = FALSE;
char  st;
char *fn,
     *old_fn;       /* the pointer fn is modified via fn++, old_fn
                       just makes it possible to free it */
FILE *infile;
int anz;
 int k, c;

file_name = B->fn;
B->fn = NULL;
	/*------------------------------------------------------------*\
	| Open input file											 |
	\*------------------------------------------------------------*/
if ( file_name == NULL )
	infile = stdin;
else
	if ( (infile = fopen (file_name, "r")) == NULL ) {
		fprintf (stderr, "get_zentr: Could not open input-file %s\n", file_name);
		exit (4);
		}
c=fscanf (infile, "%*[ \t\n\r]");
st = getc(infile);
if ( st != '#' ) {
	anz = 1;
	ungetc(st,infile);
	}
else
	c=fscanf (infile, "%u", &anz);
/*--------------------------------------------------------------------*\
|  read the matrices                                                   |
\*--------------------------------------------------------------------*/
if(anz != 0)
  B->grp->zentr = (matrix_TYP **)malloc(anz*sizeof(matrix_TYP *));
B->grp->zentr_no = anz;
for ( k = 0; k < anz; k++) {
	B->grp->zentr[k] = fget_mat(infile);
	}

	/*------------------------------------------------------------*\
	| read  file with other almost decomposable bravais-group      |
	\*------------------------------------------------------------*/
old_fn = fn = (char *) malloc(1024 *sizeof(char));
c=fscanf (infile, "%[ \t\n]", fn);
c=fscanf (infile, "%[^\n]", fn);
while(fn != NULL && fn[0] == ' ')
  fn++;

/* added free(old_fn), tilman 7/5/97 */
if(fn[0] == '\n' ||
   fn[0] == '0'  ||
   fn[0] == EOF  ||
   fn[0] == '\t' ||
   fn[0] == '%'  ||
   fn[0] == '\f' ||
   fn[0] == '\r' ||
   fn[0] == '\v'){
   free(old_fn);
   fn = NULL;
}

if(fn != NULL)
  strtok (fn, "%");
if(fn != NULL)
{                               
  B->fn = calloc( 1024, sizeof(char) );
  get_data_dir(B->fn, "tables/");
  strcat(B->fn, fn);
  free ( fn );
}
	/*------------------------------------------------------------*\
	| Close input file												|
	\*------------------------------------------------------------*/
if ( infile != stdin )
	fclose (infile);

/* inserted tilman 7/5/97 */
if (file_name != NULL) free(file_name);

}
/*{{{}}}*/
