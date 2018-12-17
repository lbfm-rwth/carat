#include"typedef.h"
#include"matrix.h"
#include"longtools.h"
#include"polyeder.h"
#include"getput.h"

/************************************************************************** \
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: refine_polyeder.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


static int g_option;

static void re_sort(F,old_no, count, test, l)
polyeder_TYP *F;
int old_no, *test, l;
int *count;
{
  int i;
  int d=0;
  for(i=0;i<old_no;i++)
  {
    if(count[i]< 0)
      free_vertex(&(F->vert[i]));
  }
  for(i=0;d<F->vert_no;i++)
  {
    while(d<F->vert_no && F->vert[d] == NULL)
      d++;
    if(d != i)
    {
      F->vert[i] = F->vert[d];
      F->vert[d] = NULL;
    }
    d++;
  }
  while(F->vert[i-1] == NULL)
    i--;
  F->vert_no = i;

  d=0;
  for(i=0;i<F->wall_no;i++)
  {
    if(test[i] == -1)
      free_wall(&(F->wall[i]));
  }
  for(i=0;d<F->wall_no;i++)
  {
    while(d<F->wall_no && F->wall[d] == NULL)
      d++;
    if(d!= i)
    {
      F->wall[i] = F->wall[d];
      F->wall[d] = NULL;
    }
    d++;
  }
  F->wall_no = i;
}


static void renumerate(v,test)
vertex_TYP *v;
int *test;
{
  int i,a;
  for(i=0;i<v->wall_no;i++)
  {
     a = v->wall[i];
     v->wall[i] = test[a];
  }
}


static void add_wall_to_vertex(i, v)
int i;
vertex_TYP *v;
{
  if(v->wall_SIZE == 0)
  {
    v->wall = (int *)malloc(extsize1 *sizeof(int));
    v->wall_SIZE = extsize1;
  }
  if(v->wall_no == v->wall_SIZE)
  {
    v->wall_SIZE += extsize1;
    v->wall = (int *)realloc(v->wall,v->wall_SIZE *sizeof(int));
  }
  v->wall[v->wall_no] = i;
  v->wall_no++;
}


static void delete_walls_from_vertex(v, test)
vertex_TYP *v;
int *test;
{
 int i,d;

 d=0;
 for(i=0; d<v->wall_no;i++)
 {
   while(d< v->wall_no && test[v->wall[d]] == -1)
     d++;
   v->wall[i] = v->wall[d];
   d++;
 }
 v->wall_no = i;
}
  


static void add_wall_to_polyeder(F, w)
polyeder_TYP *F;
wall_TYP *w;
{
  if(F->wall_SIZE == 0)
  {
    F->wall = (wall_TYP **)malloc(extsize1 *sizeof(wall_TYP *));
    F->wall_SIZE = extsize1;
  }
  if(F->wall_no == F->wall_SIZE)
  {
    F->wall_SIZE += extsize1;
    F->wall = (wall_TYP **)realloc(F->wall,F->wall_SIZE *sizeof(wall_TYP *));
  }
  F->wall[F->wall_no] = copy_wall(w);
  F->wall_no++;
}



static void add_vertex_to_polyeder(F, w)
polyeder_TYP *F;
vertex_TYP *w;
{
  if(F->vert_SIZE == 0)
  {
    F->vert = (vertex_TYP **)malloc(extsize1 *sizeof(vertex_TYP *));
    F->vert_SIZE = extsize1;
  }
  if(F->vert_no == F->vert_SIZE)
  {
    F->vert_SIZE += extsize1;
    F->vert = (vertex_TYP **)realloc(F->vert,F->vert_SIZE *sizeof(vertex_TYP *));
  }
  F->vert[F->vert_no] = w;
  F->vert_no++;
}


static vertex_TYP *gis_neighbour(i,j,F)
int i,j;
polyeder_TYP *F;
{
  int k,l,m;
  int a, u, tester1;
  int anz;
  vertex_TYP *erg;
  matrix_TYP *A;


  k=F->vert[j]->wall_no;
  if(k>F->vert[i]->wall_no)
    k=F->vert[i]->wall_no;
  erg = init_vertex(F->vert[0]->dim, k+1);
  anz=0;
  u=0;
  for(k=0;k<F->vert[i]->wall_no;k++)
  {
    a = F->vert[i]->wall[k];
    for(l=u; l<F->vert[j]->wall_no && F->vert[j]->wall[l]< a; l++);
    u=l;
    if(u<F->vert[j]->wall_no && F->vert[j]->wall[u] == a)
       { erg->wall[anz] = a; anz++;}
  }
  erg->wall_no = anz;
  if(anz < F->vert[0]->dim-2)
   { free_vertex(&erg); erg=NULL; return(NULL);}

  A = init_mat(anz,F->vert[0]->dim, "l");
  for(k=0;k<anz;k++)
    for(l=0; l<F->vert[0]->dim;l++)
      A->array.SZ[k][l] = F->wall[erg->wall[k]]->gl[l];
  a = long_row_gauss(A);
  free_mat(A);
  if(a < F->vert[0]->dim-2)
   { free_vertex(&erg); erg=NULL; return(NULL);}
  return(erg);
}



static vertex_TYP *is_neighbour(i,j,F, vertex_no)
int i,j;
polyeder_TYP *F;
int vertex_no;
{
  int k,l,m;
  int a, u, tester1;
  int anz;
  vertex_TYP *erg;


  k=F->vert[j]->wall_no;
  if(k>F->vert[i]->wall_no)
    k=F->vert[i]->wall_no;
  erg = init_vertex(F->vert[0]->dim, k+1);
  anz=0;
  u=0;
  for(k=0;k<F->vert[i]->wall_no;k++)
  {
    a = F->vert[i]->wall[k];
    for(l=u; l<F->vert[j]->wall_no && F->vert[j]->wall[l]< a; l++);
    u=l;
    if(u<F->vert[j]->wall_no && F->vert[j]->wall[u] == a)
       { erg->wall[anz] = a; anz++;}
  }
  erg->wall_no = anz;
  if(anz < F->vert[0]->dim-2)
   { free_vertex(&erg); erg=NULL; return(NULL);}

  for(k=0; k<vertex_no ;k++)
  {
   if(k != i && k != j)
   {
     tester1 = TRUE;
        u=0;
     for(l=0; l<anz && tester1 == TRUE; l++)
     {
        a = erg->wall[l];
        m=u;
        while(m<F->vert[k]->wall_no && F->vert[k]->wall[m] <a)
           m++;
        if(m == F->vert[k]->wall_no || F->vert[k]->wall[m] != a)
           tester1 = FALSE;
        u=m;
     }
     if(tester1 == TRUE)
       { free_vertex(&erg); erg=NULL; return(NULL);}
   }
  }
  return(erg);
}




/**************************************************************************\
@---------------------------------------------------------------------------
@ int refine_polyeder(F, h)
@ polyeder_TYP *F;
@ wall_TYP *h;
@ 
@  calculates the intersection of the linear Polyeder 'F' with the halfspace
@  defined by the inequality
@              H :=  h->gl * X^{tr} >= 0
@  and stores it in F.
@  If the intersection of F and H is already F, 'refine_polyeder returns 0,
@  otherwise 1.
@  All vertices of the intersecting polyeder are calculated.
@  If the intersection has not full dimension, i.e the new Polyeder is
@  non-degenerate, F->is_denerate is set to 0.
@
@  CAUTION: The result is correct only if the new Polyeder is non-degenerate.
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
int refine_polyeder(F, h)
polyeder_TYP *F;
wall_TYP *h;
{
 int i,j,k;
 int l;
 int old_no;
 int *count;
 int waste;
 int tester, *test;
 int anz;
 int p,n, z;
 vertex_TYP *v;
 int *wall_test, *necessary;

 

 p=0;n=0;
 old_no = F->vert_no;
 count = (int *)malloc(F->vert_no *sizeof(int *));
 for(i=0;i<F->vert_no;i++)
  count[i] = 0;
 for(i=0;i<old_no;i++)
 {
     count[i] = wall_times_vertex(h,F->vert[i]);
    if(count[i] > 0)
      p++;
    if(count[i] < 0)
      n++;
    if(count[i] == 0)
      z++;
 }
 if(n == 0)
 {
   free(count);
   return(0);
 }
 if(p == 0)
 {
   F->is_degenerate = TRUE;
 }


 g_option = FALSE;
 if(is_option('g') == TRUE)
    g_option = TRUE;
 wall_test = (int *)malloc(F->wall_no *sizeof(int));
 for(i=0;i<F->wall_no;i++)
     wall_test[i] = 2;
 for(i=0;i<old_no;i++)
 {
   if(count[i]>0)
   {
      for(j=0;j<F->vert[i]->wall_no;j++)
      {
        if(wall_test[F->vert[i]->wall[j]] == 2)
          wall_test[F->vert[i]->wall[j]] = 1;
        if(wall_test[F->vert[i]->wall[j]] == -1)
          wall_test[F->vert[i]->wall[j]] = 0;
      }
   }
   if(count[i] < 0)
   {
      for(j=0;j<F->vert[i]->wall_no;j++)
      {
        if(wall_test[F->vert[i]->wall[j]] == 2)
          wall_test[F->vert[i]->wall[j]] = -1;
        if(wall_test[F->vert[i]->wall[j]] == 1)
          wall_test[F->vert[i]->wall[j]] = 0;
      }
   }
 }
 necessary = (int *)malloc(F->vert_no *sizeof(int));
 for(i=0;i<old_no;i++)
   necessary[i] = FALSE;
 for(i=0;i<old_no;i++)
 {
  if(count[i] != 0)
  {
   k=0;
   for(j=0;j<F->vert[i]->wall_no && necessary[i] == FALSE;j++)
   {
     if( wall_test[F->vert[i]->wall[j]] == 0)
       k++;
     if(k>= F->vert[i]->dim-2)
     {
       if(count[i] > 0)
         necessary[i] = 1;
       if(count[i] < 0)
         necessary[i] = -1;
       if(count[i] == 0)
         necessary[i] = 0;
     }
   }
  }
 }

 if(F->vert[0]->dim == 2)
 {
   for(i=0;i<old_no;i++)
   {
     if(count[i] > 0)
       necessary[i] = 1;
     if(count[i] < 0)
       necessary[i] = -1;
     if(count[i] == 0)
       necessary[i] = 0;
   }
 }


 for(i=0;i<old_no;i++)
 {
  if(necessary[i] == -1)
  {  
    for(j=0;j<old_no;j++)
    {
     if(necessary[j] == 1)
     {
        if(g_option == TRUE)
          v = gis_neighbour(i, j, F);
        else
          v = is_neighbour(i, j, F, old_no);
        if(v != NULL)
        {
          for(k=0; k<v->dim; k++)
            v->v[k] = count[j] * F->vert[i]->v[k] - count[i] *F->vert[j]->v[k];
          normal_vertex(v);
          add_vertex_to_polyeder(F, v);
        }
     }
    }
  }
 }
  add_wall_to_polyeder(F, h);
  test = (int *)malloc((F->wall_no) *sizeof(int));
  l=0;
  for(i=0; i<F->wall_no-1;i++)
  {
    tester = 0;
    for(j=0;j<old_no && tester == 0; j++)
    {
     if(count[j] > 0)
       if(is_vertex_of_wallno(F->vert[j], i) == TRUE)
          tester = 1;
    }
    if(tester == 0)
     test[i] = -1;
    else
     {test[i] = l; l++;}
  }
  test[F->wall_no-1] = l;
  l++;
  if(test[0] == -1)
    F->is_closed = TRUE;

  for(i=0;i<old_no; i++)
  {
    if(count[i] == 0)
      add_wall_to_vertex(F->wall_no-1, F->vert[i]);
  }

  for(i=0;i<old_no;i++)
  {
     if(count[i]==0)
       delete_walls_from_vertex(F->vert[i],test);
     if(count[i]>=0)
       renumerate(F->vert[i],test);
  }
  for(i=old_no; i<F->vert_no;i++)
  {
     F->vert[i]->wall[F->vert[i]->wall_no] = F->wall_no-1;
     F->vert[i]->wall_no++;
     renumerate(F->vert[i],test);
  }

  re_sort(F, old_no, count, test,l);
  
   free(count);
   free(test);
   free(necessary);
   free(wall_test);
  return(1);
}
