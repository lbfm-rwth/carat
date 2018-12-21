
#include "typedef.h"
#include "tools.h"
#include "matrix.h"
#include "getput.h"
#include "longtools.h"
#include "datei.h"

/************************************************************************** \
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: get_symbol.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ symbol_out *get_symbol (file_name)
@ char *file_name;
@
@    This is a function for the 'Datei' programm.
@    It reads a group (bravais_TYP) from 'file_name' 
@    and additionally the name of another file stated in the
@    file 'file_name' beyond the group
@
@---------------------------------------------------------------------------
@
\**************************************************************************/

symbol_out *get_symbol (file_name)
char *file_name;

{  

char string[80];
char *str, *strin;
char merk[10];
char *fn;
symbol_out *sy;
boolean header = FALSE;
char  st;
char **piece;
FILE *infile;
int anz, teile;
 int i, j, k, l, m, c;
int gen_no, form_no, zentr_no, normal_no, cen_no;

	/*------------------------------------------------------------*\
	| Open input file											 |
	\*------------------------------------------------------------*/
if ( file_name == NULL )
	infile = stdin;
else
	if ( (infile = fopen (file_name, "r")) == NULL ) {
		fprintf (stderr, "get_symbol: Could not open input-file %s\n", file_name);
		exit (4);
		}
        gen_no = 0;
        form_no = 0;
        zentr_no = 0;
        normal_no = 0;
        cen_no = 0;

  /*--------------------------------------------------*\
  |  read header line                                  |
  \*--------------------------------------------------*/
c=fscanf (infile, "%*[ \t\n\r]");
st = getc(infile);
if ( st != '#' ) {
	gen_no = 1;
	ungetc(st,infile);
	}
else
{

c=fscanf (infile, "%[ \t\n]", string);
c=fscanf (infile, "%[^\n]",string);
strtok (string, "%");

  if ( (str = strpbrk (string, "gfznc")) == NULL )
        c=sscanf(string, "%d", &gen_no);
  else
  {
     while((str = strpbrk(str, "gfznc")) != NULL)
     {
        i = strcspn(str, "g");
        j = strcspn(str, "f");
        k = strcspn(str, "z");
        l = strcspn(str, "n");
        m = strcspn(str, "c");

        if(i< j && i<k &&i<l && i<m)
	   c=sscanf ( ++str, "%d", &gen_no);
        if(j<i && j<k && j<l && j<m)
	   c=sscanf ( ++str, "%d", &form_no);
        if(k<i && k<j && k<l && k<m)
	   c=sscanf ( ++str, "%d", &zentr_no);
        if(l<i && l<j && l<k && l<m)
	   c=sscanf ( ++str, "%d", &normal_no);
        if(m<i && m<j && m<k && m<l)
	   c=sscanf ( ++str, "%d", &cen_no);
     }
  }
}

  /*--------------------------------------------------*\
  |  read the matrices                                 |
  \*--------------------------------------------------*/
sy = (symbol_out *) malloc(sizeof(symbol_out));
sy->grp = (bravais_TYP *) malloc(sizeof(bravais_TYP));
sy->grp->gen_no = gen_no;
sy->grp->form_no = form_no;
sy->grp->zentr_no = zentr_no;
sy->grp->normal_no = normal_no;
sy->grp->cen_no = cen_no;
sy->grp->gen = (matrix_TYP **)malloc(gen_no *sizeof(matrix_TYP *));
sy->grp->form = (matrix_TYP **)malloc(form_no *sizeof(matrix_TYP *));
if (zentr_no > 0)
    sy->grp->zentr = (matrix_TYP **)malloc(zentr_no *sizeof(matrix_TYP *));
if (normal_no > 0)
    sy->grp->normal = (matrix_TYP **)malloc(normal_no *sizeof(matrix_TYP *));
if (cen_no > 0)
    sy->grp->cen = (matrix_TYP **)malloc(cen_no *sizeof(matrix_TYP *));

for ( k = 0; k < gen_no; k++)
	sy->grp->gen[k] = fget_mat(infile);
for ( k = 0; k < form_no; k++)
	sy->grp->form[k] = fget_mat(infile);
for ( k = 0; k < zentr_no; k++)
	sy->grp->zentr[k] = fget_mat(infile);
for ( k = 0; k < normal_no; k++)
	sy->grp->normal[k] = fget_mat(infile);
for ( k = 0; k < cen_no; k++)
	sy->grp->cen[k] = fget_mat(infile);
   
	/*------------------------------------------------------------*\
	| read group order                                             |
	\*------------------------------------------------------------*/
c=fscanf (infile, "%[ \t\n]", string);
c=fscanf (infile, "%[^\n]",string);
if ( *string == '%' )
  strin= NULL;
else
  strin = strtok (string, "%");
for(i=0; i<100; i++)
  sy->grp->divisors[i] = 0;
if( (strlen(strin)) != 0 && strin != NULL)
{
  i = strcspn(strin, "=");
  while(i != 0)
  {
    while(strin[0] == ' ')
       strin++;
    if((strcspn(strin, "*")) == 0)
       strin++;
    while(strin[0] == ' ')
       strin++;
    c=sscanf(strin, "%d", &j);
    itoasc(j, merk);
    k = strlen(merk);
    strin = strin+k;
    while(strin[0] == ' ')
       strin++;
    if((strcspn(strin, "^")) != 0)
      sy->grp->divisors[j] = 1;
    else
    {
      strin++;
      while(strin[0] == ' ')
         strin++;
        c=sscanf(strin, "%d", &k);
        sy->grp->divisors[j] = k;
        itoasc(k, merk);
        k = strlen(merk);
        strin = strin+k;
    }
    while(strin[0] == ' ')
       strin++;
    i = strcspn(strin, "=");
   }
   if ( (str = strpbrk (strin, "=")) != NULL )
	c=sscanf ( ++str, "%d", &sy->grp->order);
   else
     sy->grp->order = 0;
}
else
{
  sy->grp->divisors[0] = 1;
  sy->grp->order = 0;
}

	/*------------------------------------------------------------*\
	| read file-reference                                             |
	\*------------------------------------------------------------*/
sy->fn = (char *) malloc(1024 *sizeof(char));
fn = (char *) malloc(1024 *sizeof(char));
c=fscanf (infile, "%[ \t\n]", fn);
c=fscanf (infile, "%[^\n]", fn);
while(fn != NULL && fn[0] == ' ')
  fn++;
if(fn[0] == '\n' ||
    fn[0] == '0' ||
    fn[0] == EOF ||
    fn[0] == '%' ||
    fn[0] == '\f' ||
    fn[0] == '\r' ||
    fn[0] == '\v' ||
    fn[0] == '\t'){
  /* added the free: tilman 7/5/97 */
  free(fn);
  fn = NULL;
}

if(fn != NULL)
  strtok (fn, "%");
if(fn != NULL)
{
  get_data_dir(sy->fn, "tables/");
  strcat(sy->fn, fn);
}
else{
  /* added the free: tilman 7/5/97 */
  free(sy->fn);
  sy->fn = NULL;
}
	/*------------------------------------------------------------*\
	| close input file                                             |
	\*------------------------------------------------------------*/
if ( infile != stdin )
	fclose (infile);
if(sy->grp->gen_no != 0)
  sy->grp->dim = sy->grp->gen[0]->cols;

/* inserted 7/5/97 tilman (untill return (sy)) */
long_rein_formspace(sy->grp->form,sy->grp->form_no,1);
for (i=0;i<sy->grp->gen_no;i++) Check_mat(sy->grp->gen[i]);
for (i=0;i<sy->grp->form_no;i++) Check_mat(sy->grp->form[i]);
for (i=0;i<sy->grp->zentr_no;i++) Check_mat(sy->grp->zentr[i]);
for (i=0;i<sy->grp->normal_no;i++) Check_mat(sy->grp->normal[i]);
if (fn != NULL) free(fn);

return ( sy );
}
/*{{{}}}*/
