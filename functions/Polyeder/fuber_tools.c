#include "typedef.h"
#include "tools.h"
#include "matrix.h"
#include "polyeder.h"

vertex_TYP *
init_vertex_fuber (int dim, int wall_no)
{
  int j;
  vertex_TYP *erg;
  erg = (vertex_TYP *)malloc(sizeof(vertex_TYP));
  if(dim!=0)
    erg->v = (int *)malloc(dim *sizeof(int));
  erg->dim = dim;
  erg->wall_no = wall_no;
  if( wall_no > 0)
  {
    j = wall_no/extsize1 +1;
    erg->wall_SIZE = j * extsize1;
    erg->wall = (int *)malloc(erg->wall_SIZE *sizeof(int));
  }
  else
    erg->wall_SIZE=0;
    
  return(erg);
}

wall_TYP *
init_wall_fuber (int dim)
{
  wall_TYP *erg;

  erg = (wall_TYP *)malloc(sizeof(wall_TYP));
  erg->gl = (int *)malloc(dim *sizeof(int));
  erg->norm = 0;
  erg->dim = dim;
  erg->product = NULL;
  erg->nproduct = 0;
  erg->mat = NULL;
  return(erg);
}

/* 
 Aenderung: 12.10.94 JK
 Speicherplatz wieder freigeben
*/
void free_fund_domain(fund_domain* fudo)
{
int i;
for (i=0;i<fudo->wall_no;i++)
    free_wall_fuber( &fudo->wall[i]);
if (fudo->wall_no>0) 
    free( (int *)fudo->wall);
for (i=0;i<fudo->vert_no;i++)
    free_vertex_fuber( &fudo->vert[i]);
if (fudo->vert_no>0)
    free( (int *)fudo->vert);
free( (int *)fudo);
}

fund_domain *init_fund_domain(int vert_no, int wall_no)
{
  int j;
  fund_domain *erg;

  erg = (fund_domain *)malloc(sizeof(fund_domain));
  erg->vert_no = vert_no;
  erg->wall_no = wall_no;
  j = vert_no/extsize1 +1;
  j *= extsize1;
  erg->vert_SIZE = j;
  erg->vert = (vertex_TYP **)malloc(j*sizeof(vertex_TYP *));
  j = wall_no/extsize1 +1;
  j *= extsize1;
  erg->wall_SIZE = j;
  erg->wall = (wall_TYP **)malloc(j*sizeof(wall_TYP *));
  return(erg);
}




fund_domain *get_fund_domain(char *file_name)
{
int vertno, wallno;
int dim,wn,i,j;
fund_domain *F;
FILE *infile;


	/*------------------------------------------------------------*\
	| Open input file											 |
	\*------------------------------------------------------------*/
if ( file_name == NULL )
	infile = stdin;
else
	if ( (infile = fopen (file_name, "r")) == NULL ) {
		fprintf (stderr, "Could not open input-file %s\n", file_name);
		exit (4);
		}
  /*--------------------------------------------------*\
  |  read fundamental domain                                  |
  \*--------------------------------------------------*/
fscanf (infile, "%d", &vertno);
fscanf (infile, "%d", &wallno);
F = init_fund_domain(vertno, wallno);
for(i=0;i<vertno;i++)
{
  fscanf (infile, "%d", &dim);
  fscanf (infile, "%d", &wn);
  F->vert[i] = init_vertex_fuber(dim, wn);
  for(j=0;j<dim;j++)
    fscanf(infile, "%d", &F->vert[i]->v[j]); 
  for(j=0;j<wn;j++)
     fscanf(infile, "%d", &F->vert[i]->wall[j]);
}
for(i=0;i<wallno;i++)
{
  fscanf(infile, "%d", &dim);
  F->wall[i] = init_wall_fuber(dim);
  for(j=0;j<dim;j++)
    fscanf(infile, "%d", &F->wall[i]->gl[j]); 
}

   
	/*------------------------------------------------------------*\
	| close input file                                             |
	\*------------------------------------------------------------*/
if ( infile != stdin )
	fclose (infile);
return ( F );
}


void put_fund_domain(fund_domain *F)
{
  int i,j;
  printf("%d  %d\n", F->vert_no, F->wall_no);
  printf("\n");
  for(i=0;i<F->vert_no;i++)
  {
    printf("%d  %d\n", F->vert[i]->dim, F->vert[i]->wall_no);
    for(j=0;j<F->vert[i]->dim;j++)
      printf("%d  ", F->vert[i]->v[j]);
    printf("\n");
    for(j=0;j<F->vert[i]->wall_no;j++)
       printf("%d  ", F->vert[i]->wall[j]);
    printf("\n");
  }
  printf("\n");
  for(i=0;i<F->wall_no;i++)
  {
    printf("%d\n", F->wall[i]->dim);
    for(j=0;j<F->wall[i]->dim;j++)
      printf("%d  ", F->wall[i]->gl[j]);
    printf("\n");
  }
    
  fflush(stdout);
}


int 
wall_times_vertex_fuber (wall_TYP *w, vertex_TYP *v)
{
  int i;
  int e;
  e=0;
  for(i=0;i<v->dim;i++)
     e += v->v[i] * w->gl[i];
  return(e);
}

void 
free_vertex_fuber (vertex_TYP **v)
{
 if((*v)!= NULL)
 {
   free((int *)(*v)->v);
   if((*v)->wall_SIZE != 0)
     free((int *)(*v)->wall);
   free((int *)(*v));
   (*v)=NULL;
 }
}

void 
free_wall_fuber (wall_TYP **v)
{
 if((*v)!= NULL)
 {
   free((*v)->gl);
   if((*v)->mat != NULL)
      free_mat((*v)->mat);
   if ((*v)->nproduct!=0)
      free((*v)->product);
   free( (int *)(*v));
   (*v)=NULL;
 }
}

/*
 commented out this function tilman 21/3/97 
wall_TYP *mat_to_wall(M)
matrix_TYP *M;
{
  int i;
  wall_TYP *erg;
  erg = init_wall_fuber(M->cols);
  for(i=0;i<erg->dim;i++)
    erg->gl[i] = M->array.SZ[0][i];
  return(erg);
}
*/

void 
wall_standard (wall_TYP *v)
{
  int i,j;
  int w1;
 
  if(v->dim>0)
  {
    w1 = 0;
    i=0;
    while(i<v->dim && v->gl[i] == 0)
       i++;
    if(i<v->dim)
      w1 = v->gl[i];
    for(j=i+1;j<v->dim;j++)
    {
     if(v->gl[j] != 0)
       w1 = GGT(w1, v->gl[j]);
    }
    if(w1 < 0)
      w1 = -w1;
    if(w1 != 0)
    {
       for(j=i;j<v->dim;j++)
         v->gl[j] /= w1;
    }
    v->norm = 0;
    for(i=0;i<v->dim;i++)
    v->norm += v->gl[i] * v->gl[i];
  }
}
/*{{{}}}*/
