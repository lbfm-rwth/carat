#include"typedef.h"
#include"polyeder.h"
#include"matrix.h"
#include"tools.h"

/************************************************************************** \
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: polyeder_tools.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ vertex_TYP *init_vertex(dim, wall_no) 
@ int dim;
@ int wall_no;
@
@  'init_vertex' allocates a vertex_TYP *v.
@  For v the following is allocated:
@     v->dim = dim
@     v->wall_no = wall_no
@     v->v:    pointer to integer, size: dim.
@     v->wall_size: 0, if wall_no = 0
@                   (wall_no/extsize1 +1) * extsize1, else
@     v->wall: pointer to integer
@               size: 0, if wall_no = 0
@                        else (wall_no/extsize1 +1) * extsize1
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
vertex_TYP *init_vertex(dim, wall_no) 
int dim;
int wall_no;
{
  int i,j;
  vertex_TYP *erg;
  erg = (vertex_TYP *)malloc(sizeof(vertex_TYP));
  if(dim!=0)
  {
    erg->v = (int *)malloc(dim *sizeof(int));
    for(i=0; i<dim; i++)
    erg->v[i] = 0;
  }
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

/* anne, 16/07/97 */
/**************************************************************************\
@---------------------------------------------------------------------------
@ word_TYP *init_word(dim)
@ int dim;
@ 
@   'init_word' allocates a word_TYP *word.
@   For word the following is allocated:
@   word->word = NULL;
@   matrix_TYP *trans = init_mat(dim,1,"");
@---------------------------------------------------------------------------
@
\**************************************************************************/
  word_TYP *init_word(dim)
  int dim;
  {
    word_TYP *erg;
   
    erg = (word_TYP*)malloc(sizeof(word_TYP));
    erg->dim = dim;
    erg->word = NULL;
    erg->trans = NULL; 
/*  erg->trans = init_mat(dim,1,""); */
    
    return(erg);
  }

/**************************************************************************\
@---------------------------------------------------------------------------
@ wall_TYP *init_wall(dim)
@ int dim;
@   'init_wall' allocates a wall_TYP *w.
@   For w the following is allocated:
@      w->dim = dim
@      w->gl:    pointer to integer, size: dim.
@      w->mat = NULL
@      w->product = NULL
@      w->nproduct = 0
@      w->word = NULL
@      w->next_no = 0
@      w->next = NULL
@      w->ext_no = 0;
@      w->extra= NULL;
@      w->neu = 0;
@      w->paar = 0;
@---------------------------------------------------------------------------
@
\**************************************************************************/
wall_TYP *init_wall(dim)
int dim;
{
  int i,j;
  wall_TYP *erg;

  erg = (wall_TYP *)malloc(sizeof(wall_TYP));
  erg->gl = (int *)calloc(dim ,sizeof(int));
  erg->dim = dim;
  erg->mat = NULL;
  erg->product = NULL;
  erg->nproduct = 0;
  erg->word = NULL;
  erg->next_no = 0;	/* next 5 lines anne 8/10/97 */
  erg->next = NULL;
  erg->ext_no = 0;
  erg->extra = NULL;
  erg->neu = 0;
  erg->paar = 0;
  return(erg);
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ corner_TYP *init_corner()
@
@
@
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ polyeder_TYP *init_polyeder(vert_no, wall_no)
@ int vert_no, wall_no;
@ 
@   'init_polyeder' allocates a polyeder_TYP *P.
@   For P the following is allocated:
@     P->vert_no = vert_no;
@     P->wall_no = wall_no;
@     P->is_closed = FALSE;
@     P->vert_SIZE = (vert_no/extsize1 +1) * extsize1;
@     P->wall_SIZE = (wall_no/extsize1 +1) * extsize1;
@     P->corner = NULL;
@     P->corner_no = 0, P->corner_SIZE = 0;
@     P->vert = **vertex_TYP, size:  P->vert_no.
@     P->wall = **wall_TYP,   size:  P->wall_no.
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
polyeder_TYP *init_polyeder(vert_no, wall_no)
int vert_no, wall_no;
{
  int i,j;
  polyeder_TYP *erg;

  erg = (polyeder_TYP *)malloc(sizeof(polyeder_TYP));
  erg->vert_no = vert_no;
  erg->wall_no = wall_no;
  erg->is_closed = FALSE;
  erg->is_degenerate = FALSE;
  j = vert_no/extsize1 +1;
  j *= extsize1;
  erg->vert_SIZE = j;
  erg->vert = (vertex_TYP **)malloc(j*sizeof(vertex_TYP *));
  for(i=0;i<erg->vert_SIZE;i++)
    erg->vert[i] = NULL;
  j = wall_no/extsize1 +1;
  j *= extsize1;
  erg->wall_SIZE = j;
  erg->wall = (wall_TYP **)malloc(j*sizeof(wall_TYP *));
  for(i=0;i<erg->wall_SIZE;i++)
    erg->wall[i] = NULL;
  erg->corner = NULL;			/* next 3 lines inserted by anne */
  erg->corner_no = 0;
  erg->corner_SIZE = 0;
  return(erg);
}





/**************************************************************************\
@---------------------------------------------------------------------------
@ polyeder_TYP *get_polyeder(file_name)
@ char *file_name;
@ 
@    Read a polyeder from the file 'file_name'.
@    If 'file_name' is NULL, get+polyeder reads from standard input.
@ 
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
polyeder_TYP *get_polyeder(file_name)
char *file_name;
{  
int vertno, wallno;
int dim,wn,i,j,c;
polyeder_TYP *F;
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
c=fscanf (infile, "%d", &vertno);
c=fscanf (infile, "%d", &wallno);
F = init_polyeder(vertno, wallno);
for(i=0;i<vertno;i++)
{
  c=fscanf (infile, "%d", &dim);
  c=fscanf (infile, "%d", &wn);
  F->vert[i] = init_vertex(dim, wn);
  for(j=0;j<dim;j++)
    c=fscanf(infile, "%d", &F->vert[i]->v[j]); 
  for(j=0;j<wn;j++)
     c=fscanf(infile, "%d", &F->vert[i]->wall[j]);
}
for(i=0;i<wallno;i++)
{
  c=fscanf(infile, "%d", &dim);
  F->wall[i] = init_wall(dim);
  for(j=0;j<dim;j++)
    c=fscanf(infile, "%d", &F->wall[i]->gl[j]); 
}
  c=fscanf(infile, "%d", &F->is_closed);
  c=fscanf(infile, "%d", &F->is_degenerate);

   
	/*------------------------------------------------------------*\
	| close input file                                             |
	\*------------------------------------------------------------*/
if ( infile != stdin )
	fclose (infile);
return ( F );
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ void put_polyeder(F)
@ polyeder_TYP *F;
@ 
@    prints a polyeder_TYP to standard output
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
void put_polyeder(F)
polyeder_TYP *F;
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
    
  printf("\n");
  printf("%d  %d\n", F->is_closed, F->is_degenerate);
  fflush(stdout);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ int wall_times_vertex(w, v)
@ wall_TYP *w;
@ vertex_TYP *v;
@ 
@     calculates the sum of v->v[i] * w->gl[i].
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
int wall_times_vertex(w, v)
wall_TYP *w;
vertex_TYP *v;
{
  int i,e;
  e=0;
  for(i=0;i<v->dim;i++)
      e += (v->v[i] * w->gl[i]);
  return(e);
}



/**************************************************************************\
@---------------------------------------------------------------------------
@ void free_vertex(v)
@ vertex_TYP **v;
@
@    frees, what was allocated in *v and sets (*v) = NULL 
@---------------------------------------------------------------------------
@
\**************************************************************************/
void free_vertex(v)
vertex_TYP **v;
{
 int i;
 if((*v)!= NULL)
 {
   free((*v)->v);
   if((*v)->wall_SIZE != 0)
     free((*v)->wall);
   free((*v));
   (*v)=NULL;
 }
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ void free_wall(v)
@ wall_TYP **v;
@ 
@    frees, what was allocated in *v and sets (*v) = NULL 
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
void free_wall(v)
wall_TYP **v;
{
 int i;

 if((*v)!= NULL)
 {
   free((*v)->gl);
   if((*v)->mat != 0)
     free_mat((*v)->mat);
   if((*v)->next != NULL){		/* next 7 lines anne, 8/10/97 */
     for(i=0; i<(*v)->next_no; i++)
         free((*v)->next[i]);
     free((*v)->next);
   }
   if((*v)->extra != NULL){
     for(i=0; i<(*v)->ext_no; i++)
         free((*v)->extra[i]);
     free((*v)->extra);
   } 
   if((*v)->word != NULL){		/* next 7 lines anne, 1.4. 98 */
     if((*v)->word->trans != NULL){
        free_mat((*v)->word->trans); (*v)->word->trans = NULL; } 
     if((*v)->word->word != NULL)
        free((*v)->word->word);
   free((*v)->word);
   }
   free((*v));
   if((*v)->product != NULL)
     free((*v)->product);

   (*v)=NULL;
 }
}

/* anne, 16/07/97  */
/**************************************************************************\
@---------------------------------------------------------------------------
@ void free_word(word)
@ word_TYP *word; 
@
@    frees, what was allocated in word and sets word = NULL 
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
void free_word(word)
word_TYP *word;
{
   if (word != NULL){
      if (word->trans != NULL)
         free_mat(word->trans);
      if (word->word[0] > 0 &&
          word->word != NULL)
         free(word->word);
      free(word);
      word = NULL;
   }
}

/**************************************************************************\
@---------------------------------------------------------------------------
@ wall_TYP *mat_to_wall(M)
@ matrix_TYP *M;
@ 
@     creates a wall_TYP *w with w->dim = M->cols
@     and the entries of w->gl equal to the entries of the first row of M.
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
wall_TYP *mat_to_wall(M)
matrix_TYP *M;
{
  int i;
  wall_TYP *erg;
  erg = init_wall(M->cols);
  for(i=0;i<erg->dim;i++)
    erg->gl[i] = M->array.SZ[0][i];
  return(erg);
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ void normal_wall(v)
@ wall_TYP *v;
@ 
@     Divides the entries of v->gl by their greatest common divisor.
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
void normal_wall(v)
wall_TYP *v;
{
  int i,j;
  int w1,w2;
 
  if(v->dim>0)
  {
    i=0;
    while(i<v->dim && v->gl[i] == 0)
       i++;
    if(i<v->dim)
      w1 = v->gl[i];
    for(j=i+1;j<v->dim && w1 != 1 && w1 != -1;j++)
    {
     if(v->gl[j] != 0)
     {
       w2 = GGT(w1, v->gl[j]);
       w1 = w2;
     }
    }
    if(w1 < 0)
      w1 = -w1;
    if(w1 != 0)
    {
       for(j=0;j<v->dim;j++)
         v->gl[j] /= w1;
    }
  }
}




/**************************************************************************\
@---------------------------------------------------------------------------
@ void normal_vertex(v)
@ vertex_TYP *v;
@ 
@     Divides the entries of v->v by their greatest common divisor.
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
void normal_vertex(v)
vertex_TYP *v;
{
  int i,j;
  int w1,w2;
 
  if(v->dim>0)
  {
    i=0;
    while(i<v->dim && v->v[i] == 0)
       i++;
    if(i<v->dim)
      w1 = v->v[i];
    for(j=i+1;j<v->dim && w1 != 1 && w1 != -1;j++)
    {
     if(v->v[j] != 0)
     {
       w2 = GGT(w1, v->v[j]);
       w1 = w2;
     }
    }
    if(w1 < 0)
      w1 = -w1;
    if(w1 != 0)
    {
       for(j=0;j<v->dim;j++)
         v->v[j] /= w1;
    }
  }
}





/**************************************************************************\
@---------------------------------------------------------------------------
@ int is_vertex_of_wallno(v, w)
@ vertex_TYP *v;
@ int w;
@ 
@     Checks if the inter w is an entry of v->wall.
@     If w is an entry of w->wall, the result is 1, otherwise 0.
@     CAUTION: the entries of v->wall have to be ordered, t.m. 
@              v->wall[i] < v->wall[i+1].
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
int is_vertex_of_wallno(v, w)
vertex_TYP *v;
int w;
{
  int o,u,t;
  o=w;
  u=0;
  if(o>= v->wall_no)
    o=v->wall_no-1;
  while(o>u)
  {
    t=(o+u)/2;
    if(v->wall[t] < w)
      u=t+1;
    else
      o=t;
  } 
  if(w == v->wall[u])
    return(TRUE);
  return(FALSE);
}

/**************************************************************************\
@---------------------------------------------------------------------------
@ word_TYP *copy_word(w)
@ word_TYP *w;
@  make returns a copy of w.
@ in w->word[0] the length of w->word is encoded.
@---------------------------------------------------------------------------
\**************************************************************************/
word_TYP *copy_word(w)
word_TYP *w;
{
  word_TYP *erg;
  int	   i;

  erg = init_word(w->dim);
  erg->word = (int*)malloc((w->word[0]+1)* sizeof(int));
  for(i=0; i<= w->word[0]; i++)
     erg->word[i] = w->word[i];
  if(w->trans != NULL)
     erg->trans = copy_mat(w->trans);

  return(erg);
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ wall_TYP *copy_wall(w)
@ wall_TYP *w;
@  make returns a copy of w.
@---------------------------------------------------------------------------
@
\**************************************************************************/
wall_TYP *copy_wall(w)
wall_TYP *w;
{
  wall_TYP *erg;
  int i;
  int j;

  erg = init_wall(w->dim);
  for(i=0;i<w->dim;i++)
   erg->gl[i] = w->gl[i];
  if(w->mat != NULL)
   erg->mat = copy_mat(w->mat);
  else
   erg->mat = NULL;
  if(w->word != NULL)                    /*anne, 16/07/97 */
   erg->word = copy_word(w->word);
  else
   erg->word = NULL;
  if(w->next != NULL){			/* 3 lines anne 8/10/97 */
   erg->next_no = w->next_no;
   erg->next = (int **)malloc(w->next_no *sizeof(int*));
   for(j=0; j<w->next_no; j++){
      erg->next[j] = (int*)malloc(w->dim * sizeof(int));
      for(i=0;i<w->dim;i++)
         erg->next[j][i] = w->next[j][i]; 
   }
  }
  else{
   erg->next = NULL;
   erg->next_no = 0;
  }
  if(w->extra != NULL){		/* 12 lines, anne 14/10/97 */
   erg->ext_no = w->ext_no;
   erg->extra = (int **)malloc(w->ext_no *sizeof(int*));
   for(j=0; j<w->ext_no; j++){
      erg->extra[j] = (int*)malloc(w->dim * sizeof(int));
      for(i=0;i<w->dim;i++)
         erg->extra[j][i] = w->extra[j][i]; 
   }
  }
  else{
   erg->extra = NULL;
   erg->ext_no = 0;
  }
  erg->neu = w->neu;
  erg->nproduct = w->nproduct;
  if(w->nproduct != 0)
  {
     if( (erg->product = (int *)malloc(w->nproduct *sizeof(int))) == 0)
     {
        printf("malloc failed in copy_wall\n");
        exit(2);
     }
  }
  for(i=0;i<w->nproduct;i++)
    erg->product[i] = w->product[i];
  return(erg);
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ void free_polyeder(P)
@ polyeder_TYP *P;
@
@  frees all the pointers allocated in P, and P itself
@---------------------------------------------------------------------------
@
\**************************************************************************/
void free_polyeder(P)
polyeder_TYP *P;
{
  int i;
  for(i=0;i<P->vert_no;i++)
  {
    if(P->vert[i] != NULL)
      free_vertex(&P->vert[i]);
  }
  free(P->vert);
  for(i=0;i<P->wall_no; i++)
  {
   if(P->wall[i] != 0)
        free_wall(&P->wall[i]);
  }
  if(P->corner != NULL) 
    free(P->corner);
  if(P->wall != NULL)
    free(P->wall);
  free(P);
}
