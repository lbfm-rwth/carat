#include "typedef.h"    

#include "matrix.h"
#include "getput.h"
#include "datei.h"
#include "longtools.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: read_symbol.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
static void clear(merk)
char *merk;
{
  int i, k;
  k = strlen(merk);
  for(i=0; i<k; i++)
    merk[i] = '\0';
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ symbol_out *read_symbol(file_name)
@ char *file_name;
@
@   reads the symbol for the Crystal-fammily used in the 'Datei' programm
@   and constructs the first almost decomposable group by reading the
@   atoms from the libary.
@---------------------------------------------------------------------------
@
\**************************************************************************/
symbol_out *read_symbol(file_name)
char *file_name;

{
char  string[80],
      slash, *str ;
char f[1024], dat[1024];
/* changed from "static char fn[80]" to switch away from static variables */
char *fn;
FILE *infile;
int i, j, k, l, m, n, p, q, x, groesser;
int no, c;
int breite;
char merk[15];
char merk1[10];
int index;
char konst[MAXDIM][80];
int konst_dim;
int komp[MAXDIM];
int artgleich[MAXDIM];
int dim = 0, konstit = 0;
int zerleg[MAXDIM][5];
bravais_TYP **grps;
symbol_out *erg;
matrix_TYP *In;

/* it is assumed that string contains a 0-string, all not all
   compiler will asure that (inserted 24/4/97 tilman) */
for (i=0;i<80;i++) string[i] = 0;

/* inserted to switch away from static variables */
fn = (char *) calloc(1024,sizeof(char));

/*------------------------------------------------------------*\
| Open input file               	         	       |
\*------------------------------------------------------------*/
if ( file_name == NULL )
	infile = stdin;
else
	if ( (infile = fopen (file_name, "r")) == NULL ) {
		fprintf (stderr, "read_symbol: Could not open input-file %s\n", file_name);
		exit (4);
		}
/*------------------------------------------------------------*\
| Read and scan header line										|
\*------------------------------------------------------------*/
for( i=0; i<MAXDIM; i++)
  komp[i] = 0;
for(i=0; i<MAXDIM; i++)
  for(j=0; j<5;j++)
    zerleg[i][j] = 0;
printf("Please input the symbol for the crystal-family: ");
c=fscanf (infile, "%[ \t\n]", string);
c=fscanf (infile, "%[^\n]",string);
strtok (string, "%");

right_order(string);

printf("%s\n",string);

str = string;
while( strlen(str) != 0)
{
  i = strcspn(str, ";,");
  if(i>0)
  {
    sscanf(str, "%d", &konst_dim);
    dim += konst_dim;
    if(dim > MAXDIM)
    {
       printf("dimension of a.d. bravais-group is to big\n");
       exit(3);
    }
    zerleg[konstit][0] = konst_dim;
    itoasc(konst_dim, konst[konstit]);
    k = strcspn(str, ";");
    if(k == i)
      zerleg[konstit][3] = 1;
    j = strcspn(str, "-");
    l = strcspn(str, "'");
    if(j<i)
    {
       strcat(konst[konstit], "-");
       str = str+j+1;
       sscanf(str, "%d", &index);
       zerleg[konstit][1] = index;
       itoasc(index, merk);
       strcat(konst[konstit], merk);
       clear(merk);
       if(l<i)
       {
         strcat(konst[konstit], "'");
         zerleg[konstit][2] = 1;
       }
       i = i-j-1;
    }
    konstit++;
  }
  str = str+i+1;
}

/*------------------------------------------------------------*\
| Close input file												|
\*------------------------------------------------------------*/
if ( infile != stdin )
	fclose (infile);

erg = (symbol_out *) malloc(sizeof(symbol_out));
erg->grp = (bravais_TYP *) calloc(1, sizeof(bravais_TYP));
erg->grp->dim = dim;

/*--------------------------------------------------------------------*\
|  swap the atoms into the right order                                 |
\*--------------------------------------------------------------------*/
for(i=0; i<MAXDIM-1; i++)
{
  if(zerleg[i][0] != zerleg[i+1][0] || zerleg[i][1] != zerleg[i+1][1] || zerleg[i][2]  != zerleg[i+1][2])
    zerleg[i][3] = 1;
}
if(zerleg[MAXDIM-1][0] != 0)
  zerleg[MAXDIM-1][3] = 1;
i=0;
while(i<MAXDIM)
{
  while( i< MAXDIM && zerleg[i][0] == 0)
    i++;
  k = 1;
  while(i<MAXDIM && zerleg[i][3] == 0 && zerleg[i][0] != 0)
  {  i++; k++;}
  if(i != MAXDIM)
    zerleg[i][4] = k;
  i++;
}
for(i=0; i<MAXDIM; i++)
{
  for(j=i+1; j<MAXDIM; j++)
  {
    groesser = 1;
    if(zerleg[i][4] == 0 && zerleg[j][4] != 0)
       groesser = 0;
    for(k=0; k<3 && zerleg[j][k] == zerleg[i][k]; k++);
    if(k == 0 && zerleg[i][k] < zerleg[j][k])
      groesser = 0;
    if(k == 1 && zerleg[i][1] > zerleg[j][1])
      groesser = 0;
    if(k == 2 && zerleg[i][2] < zerleg[j][2])
      groesser = 0;
    if(k == 3 && zerleg[i][4] < zerleg[j][4])
      groesser = 0;
    if(groesser == 0)
    {
      for(k=0; k<5; k++)
      {
         x = zerleg[i][k]; zerleg[i][k] = zerleg[j][k]; zerleg[j][k] = x;
      }
    }
  } 
}
/************************
for(i=0; i<konstit; i++)
{
  for(j=0; j<5; j++)
     printf("%d ", zerleg[i][j]);
  printf("\n");
}
*************************/

konstit = 0;
for(i=0; i<MAXDIM; i++)
{
  if(zerleg[i][4] != 0)
  {
    itoasc(zerleg[i][0], konst[konstit]);
    if(zerleg[i][1] != 0)
    {
       strcat(konst[konstit], "-");
       itoasc(zerleg[i][1], merk);
       strcat(konst[konstit], merk);
       clear(merk);
       if(zerleg[i][2] != 0)
         strcat(konst[konstit], "\'");
    }
    konstit++;
  }
}

/*--------------------------------------------------------------------*\
|  read the atoms                                                      |
\*--------------------------------------------------------------------*/
grps = (bravais_TYP **) malloc(konstit *sizeof(bravais_TYP *));
get_data_dir(dat, "tables/atoms/");
for(i=0; i<konstit; i++)
{
   strcpy(f, dat);
   itoasc(zerleg[i][0], merk);
   strcat(f, merk);
   if(zerleg[i][1] != 0)
   {
     strcat(f, "-");
     itoasc(zerleg[i][1], merk);
     strcat(f, merk);
     clear(merk);
   }
   if(zerleg[i][2] != 0)
   {
     strcat(f, "\'");
   }    
   /* fprintf(stderr,"file-name: %s\n",f); */
   grps[i] = get_bravais(f);
   for(j=0; j<grps[i]->form_no; j++)
     Check_mat(grps[i]->form[j]);
}

/*--------------------------------------------------------------------*\
|  calculate the generators of erg->grp                                |
\*--------------------------------------------------------------------*/
erg->grp->gen_no = 0;
erg->grp->form_no = 0;
erg->grp->zentr_no = 0;
erg->grp->normal_no = 0;
for(i=0; i<konstit; i++)
  erg->grp->gen_no += grps[i]->gen_no;
erg->grp->gen = (matrix_TYP **)malloc(erg->grp->gen_no *sizeof(matrix_TYP));
for(i=0; i<erg->grp->gen_no; i++)
  erg->grp->gen[i]  = init_mat(erg->grp->dim, erg->grp->dim, "");
j=0;
p=0;
for(i=0; i<konstit; i++)
{
 for(k=0; k<grps[i]->gen_no; k++)
 {
  for(l=0; l<zerleg[i][4] * grps[i]->dim; l += grps[i]->dim)
  {
   for(m=0; m<grps[i]->dim; m++)
    for(n=0; n<grps[i]->dim; n++)
     erg->grp->gen[p+k]->array.SZ[j+l+m][j+l+n]=grps[i]->gen[k]->array.SZ[m][n];
  }
  for(m=0; m<j; m++)
     erg->grp->gen[p+k]->array.SZ[m][m] = 1;
  for(m=( j + zerleg[i][4] * grps[i]->dim); m<erg->grp->dim; m++)
     erg->grp->gen[p+k]->array.SZ[m][m] = 1;
 }
   j = j + zerleg[i][4] * grps[i]->dim; 
   p = p + grps[i]->gen_no;
}

/*--------------------------------------------------------------------*\
|  calculate the invariant forms of erg->grp                           |
\*--------------------------------------------------------------------*/
for(i=0; i<konstit; i++)
{
  for(j=0; j<grps[i]->form_no; j++)
  {
     if(grps[i]->form[j]->flags.Symmetric == TRUE)
        erg->grp->form_no += (zerleg[i][4] * (zerleg[i][4] +1)/2);
     if(grps[i]->form[j]->flags.Symmetric == FALSE)
        erg->grp->form_no += (zerleg[i][4] * (zerleg[i][4] -1)/2);
  }
}
erg->grp->form = (matrix_TYP **)malloc(erg->grp->form_no *sizeof(matrix_TYP *));
for(i=0; i<erg->grp->form_no; i++)
  erg->grp->form[i]  = init_mat(erg->grp->dim, erg->grp->dim, "");
j=0;
q=0;
for(i=0; i<konstit; i++)
{
 for(k=0; k<grps[i]->form_no; k++)
 {
  for(l=j; l< j+zerleg[i][4] * grps[i]->dim; l += grps[i]->dim)
  {
   if(grps[i]->form[k]->flags.Symmetric == TRUE)
   {
    for(n=0; n<grps[i]->dim; n++) 
    {
     for(p=0; p<grps[i]->dim; p++)
     {
       erg->grp->form[q]->array.SZ[l+n][l+p] = grps[i]->form[k]->array.SZ[n][p];
     }
    }
    q++;
   }
   for(m=(l+grps[i]->dim); m<(j+zerleg[i][4] * grps[i]->dim); m +=grps[i]->dim)
   {
    for(n=0; n<grps[i]->dim; n++) 
    {
     for(p=0; p<grps[i]->dim; p++)
     {
       erg->grp->form[q]->array.SZ[l+n][m+p] = grps[i]->form[k]->array.SZ[n][p];
       erg->grp->form[q]->array.SZ[m+p][l+n] = grps[i]->form[k]->array.SZ[n][p];
     }
    }
    q++;
   }
  }
 }
 j = j + zerleg[i][4] * grps[i]->dim; 
}
/* inserted 06/05/97 tilman: assure that the full integral lattice is passed */
long_rein_formspace(erg->grp->form,erg->grp->form_no,1);

/*--------------------------------------------------------------------*\
| calculate the centralizer of the group                               |
\*--------------------------------------------------------------------*/
for(i=0; i<konstit; i++)
{
  erg->grp->cen_no = erg->grp->cen_no + grps[i]->cen_no;
  if(zerleg[i][4] >= 2)
    erg->grp->cen_no = erg->grp->cen_no + grps[i]->zentr_no + 1;
  if(zerleg[i][4] > 2)
    erg->grp->cen_no ++;
}
erg->grp->cen = (matrix_TYP **)malloc(erg->grp->cen_no *sizeof(matrix_TYP));
for(i=0; i<erg->grp->cen_no; i++)
  erg->grp->cen[i]  = init_mat(erg->grp->dim, erg->grp->dim, "");
j=0;
no = 0;
for(i=0; i<konstit; i++)
{
   for(k=0; k<grps[i]->cen_no; k++)
   {
     for(l=0; l<j; l++)
       erg->grp->cen[no]->array.SZ[l][l] = 1;
     for(l=j+grps[i]->dim; l<erg->grp->dim; l++)
       erg->grp->cen[no]->array.SZ[l][l] = 1;
     for(l=0; l<grps[i]->dim; l++)
       for(m=0; m<grps[i]->dim; m++)
        erg->grp->cen[no]->array.SZ[j+l][j+m] = grps[i]->cen[k]->array.SZ[l][m];
     no++;
   }
   if(zerleg[i][4] >=2)
   {
     for(k=0; k<grps[i]->zentr_no; k++)
     {
       for(l=0; l<j; l++)
         erg->grp->cen[no]->array.SZ[l][l] = 1;
       for(l=j; l<j+grps[i]->dim; l++)
         erg->grp->cen[no]->array.SZ[l][l] = -1;
       for(l=j+grps[i]->dim; l<erg->grp->dim; l++)
         erg->grp->cen[no]->array.SZ[l][l] = 1;
       for(l=0; l<grps[i]->dim; l++)
          for(m=0; m<grps[i]->dim; m++)
            erg->grp->cen[no]->array.SZ[j+l+grps[i]->dim][j+m] = grps[i]->zentr[k]->array.SZ[l][m];
       no++;
     }
       for(l=0; l<j; l++)
         erg->grp->cen[no]->array.SZ[l][l] = 1;
       for(l=j+(2 * grps[i]->dim); l<erg->grp->dim; l++)
         erg->grp->cen[no]->array.SZ[l][l] = 1;
       for(l=0; l<grps[i]->dim; l++)
       {
          erg->grp->cen[no]->array.SZ[j+l+grps[i]->dim][j+l] = 1;
          erg->grp->cen[no]->array.SZ[j+l][j+l+grps[i]->dim] = 1;
       }
       no++;
   }
   if(zerleg[i][4] > 2)
   {
       for(l=0; l<j; l++)
         erg->grp->cen[no]->array.SZ[l][l] = 1;
       for(l=j+(zerleg[i][4] * grps[i]->dim); l<erg->grp->dim; l++)
         erg->grp->cen[no]->array.SZ[l][l] = 1;
       for(l=0; l< ((zerleg[i][4] -1) * grps[i]->dim); l++)
          erg->grp->cen[no]->array.SZ[l+j+grps[i]->dim][l+j] = 1;
       for(l=0; l<grps[i]->dim; l++)
          erg->grp->cen[no]->array.SZ[j+l][l+j+((zerleg[i][4]-1) * grps[i]->dim)] = 1;
     no++;
   }
   
   j = j + zerleg[i][4] * grps[i]->dim; 
}

/*--------------------------------------------------------------------*\
| calculate the normalizer of the group                               |
\*--------------------------------------------------------------------*/

for(i=0; i<konstit; i++)
{
        artgleich[i] = 0;
        erg->grp->normal_no += grps[i]->normal_no;
        if(i+1 < konstit)
        {
          if(zerleg[i][0] == zerleg[i+1][0] && zerleg[i][1] == zerleg[i+1][1] &&
             zerleg[i][2] == zerleg[i+1][2] && zerleg[i][3] == zerleg[i+1][3] &&
             zerleg[i][4] == zerleg[i+1][4])
          {
             artgleich[i] = TRUE;
             erg->grp->normal_no++;
          }
        }
}
erg->grp->normal = (matrix_TYP **)malloc(erg->grp->normal_no *sizeof(matrix_TYP));
for(i=0; i<erg->grp->normal_no; i++)
  erg->grp->normal[i]  = init_mat(erg->grp->dim, erg->grp->dim, "");
j=0;
no = 0;
for(i=0; i<konstit; i++)
{
   for(k=0; k<grps[i]->normal_no; k++)
   {
      for(l=0; l<j;l++)
        erg->grp->normal[no]->array.SZ[l][l] = 1;
      for(l=(j+zerleg[i][4] * grps[i]->dim); l<erg->grp->dim; l++)
        erg->grp->normal[no]->array.SZ[l][l] = 1;
      for(l=j; l< (j+zerleg[i][4] * grps[i]->dim); l += grps[i]->dim)
      {
         for(m=0; m<grps[i]->dim; m++)
           for(p=0; p<grps[i]->dim; p++)
         erg->grp->normal[no]->array.SZ[m+l][p+l] = grps[i]->normal[k]->array.SZ[m][p];
      }
      no++;
   }
   if(artgleich[i] == TRUE)
   {
     breite = grps[i]->dim * zerleg[i][4];
     for(l=0; l<j; l++)
       erg->grp->normal[no]->array.SZ[l][l] = 1;
     for(l=(j + 2 * breite); l<erg->grp->dim; l++)
       erg->grp->normal[no]->array.SZ[l][l] = 1;
     for(l=j; l< (j+breite); l++)
     {
       erg->grp->normal[no]->array.SZ[l][l+breite] = 1;
       erg->grp->normal[no]->array.SZ[l+breite][l] = 1;
     }
     no++;
   }
   j = j + zerleg[i][4] * grps[i]->dim; 
}

/*--------------------------------------------------------------------*\
| if no additional normalizer or centraliser nessecarry, put identity  |
\*--------------------------------------------------------------------*/
if(erg->grp->normal_no == 0 || erg->grp->cen_no == 0)
{
  if(erg->grp->normal_no == 0)
  {
     erg->grp->normal_no = 1;
     erg->grp->normal = (matrix_TYP **) malloc(1 *sizeof(matrix_TYP *));
     erg->grp->normal[0] = einheitsmatrix(erg->grp->dim);
  }
  if(erg->grp->cen_no == 0)
  {
     erg->grp->cen_no = 1;
     erg->grp->cen = (matrix_TYP **) malloc(1 *sizeof(matrix_TYP *));
     erg->grp->cen[0] = einheitsmatrix(erg->grp->dim);
  }
}

/*--------------------------------------------------------------------*\
| calculate the order of erg->grp                                      |
\*--------------------------------------------------------------------*/
for(j=0; j<100; j++)
  erg->grp->divisors[j] = 0;
erg->grp->order = 1;
for(i=0; i<konstit; i++)
{
  for(j=0; j<100; j++)
    erg->grp->divisors[j] += grps[i]->divisors[j];
  erg->grp->order *= grps[i]->order;
}

/*--------------------------------------------------------------------*\
| find file where erg->grp->zentr are stored                           |
\*--------------------------------------------------------------------*/
get_data_dir(fn, "tables/dim");
itoasc(erg->grp->dim, merk);
strcat(fn, merk);
strcat(fn, "/");
for(i=0; i<konstit; i++)
{
  itoasc(zerleg[i][0], merk);
  if(zerleg[i][1] != 0)
  {
    strcat(merk, "-");
    itoasc(zerleg[i][1], merk1);
    strcat(merk, merk1);
  }
  if(zerleg[i][2] != 0)
    strcat(merk, "\'");
  for(j=0; j<zerleg[i][4]; j++)
  {
    strcat(fn, merk);
    if(j != (zerleg[i][4] -1))
      strcat(fn, ",");
    else
    {
      if(i!= (konstit -1))
        strcat(fn, ";");
    }
  }
}
erg->fn = fn;

/* printf("%s\n", erg->fn); */

 return(erg);
}
/*{{{}}}*/
