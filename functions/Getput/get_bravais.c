#include "typedef.h"
#include "getput.h"
#include "matrix.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  get_bravais.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/




/**************************************************************************\
@---------------------------------------------------------------------------
@ bravais_TYP *get_bravais (file_name)
@ char *file_name;
@    reads a bravais_TYP from the file 'file_name'
@---------------------------------------------------------------------------
@
\**************************************************************************/
bravais_TYP *get_bravais (file_name)
char *file_name;

{  

char string[1024];
char *str, *strin;
char merk[256];
bravais_TYP *grp;
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
		fprintf (stderr, "get_bravais: Could not open input-file %s\n", file_name);
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
if (c=fscanf (infile, "%[^\n]",string) == EOF) {
	*string = '\0';
}
strtok (string, "%");

  if ( (str = strpbrk (string, "gfznc")) == NULL )
        sscanf(string, "%d", &gen_no);
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
	   sscanf ( ++str, "%d", &gen_no);
        if(j<i && j<k && j<l && j<m)
	   sscanf ( ++str, "%d", &form_no);
        if(k<i && k<j && k<l && k<m)
	   sscanf ( ++str, "%d", &zentr_no);
        if(l<i && l<j && l<k && l<m)
	   sscanf ( ++str, "%d", &normal_no);
        if(m<i && m<j && m<k && m<l)
	   sscanf ( ++str, "%d", &cen_no);
     }
  }
}

  /*--------------------------------------------------*\
  |  read the matrices                                 |
  \*--------------------------------------------------*/
grp = (bravais_TYP *) malloc(sizeof(bravais_TYP));
grp->gen_no = gen_no;
grp->form_no = form_no;
grp->zentr_no = zentr_no;
grp->normal_no = normal_no;
grp->cen_no = cen_no;
if ( gen_no >0 ) {
  grp->gen = (matrix_TYP **)malloc(gen_no *sizeof(matrix_TYP *));
} else {
  grp->gen = NULL;
} 
if ( form_no>0  ) {
  grp->form = (matrix_TYP **)malloc(form_no *sizeof(matrix_TYP *));
} else {
  grp->form = NULL;
}
if ( zentr_no >0 ) {
  grp->zentr = (matrix_TYP **)malloc(zentr_no *sizeof(matrix_TYP *));
} else {
  grp->zentr = NULL;
} 
if ( normal_no>0  ) {
  grp->normal = (matrix_TYP **)malloc(normal_no *sizeof(matrix_TYP *));
} else {
  grp->normal = NULL;
}
if ( cen_no>0  ) {
  grp->cen = (matrix_TYP **)malloc(cen_no *sizeof(matrix_TYP *));
} else {
  grp->cen = NULL;
}

for ( k = 0; k < gen_no; k++)
	grp->gen[k] = fget_mat(infile);
for ( k = 0; k < form_no; k++)
	grp->form[k] = fget_mat(infile);
for ( k = 0; k < zentr_no; k++)
	grp->zentr[k] = fget_mat(infile);
for ( k = 0; k < normal_no; k++)
	grp->normal[k] = fget_mat(infile);
for ( k = 0; k < cen_no; k++)
	grp->cen[k] = fget_mat(infile);
   
	/*------------------------------------------------------------*\
	| read group order                                             |
	\*------------------------------------------------------------*/
c=fscanf (infile, "%[ \t\n]", string);
if (c=fscanf (infile, "%[^\n]",string) == EOF) {
	*string = '\0';
}
if ( *string == '%' )
  strin= NULL;
else
  strin= strtok (string, "%");
for(i=0; i<100; i++)
  grp->divisors[i] = 0;
/* changed 15/5/97 tilman from 
if( strin != NULL && (strlen(strin)) != 0 ) to */
if( strin != NULL && (strlen(strin)) != 1 )
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
    sscanf(strin, "%d", &j);
    itoasc(j, merk);
    k = strlen(merk);
    strin = strin+k;
    while(strin[0] == ' ')
       strin++;
    if((strcspn(strin, "^")) != 0)
      grp->divisors[j] = 1;
    else
    {
      strin++;
      while(strin[0] == ' ')
         strin++;
        sscanf(strin, "%d", &k);
        grp->divisors[j] = k;
        itoasc(k, merk);
        k = strlen(merk);
        strin = strin+k;
    }
    while(strin[0] == ' ')
       strin++;
    i = strcspn(strin, "=");
   }
   if ( (str = strpbrk (strin, "=")) != NULL )
	sscanf ( ++str, "%d", &grp->order);
   else
     grp->order = 0;
}
else
{
  grp->divisors[0] = 1;
  grp->order = 0;
}
   
	/*------------------------------------------------------------*\
	| close input file                                             |
	\*------------------------------------------------------------*/
if ( infile != stdin )
	fclose (infile);
if(grp->gen_no != 0)
  grp->dim = grp->gen[0]->cols;
return ( grp );
}
/*{{{}}}*/
