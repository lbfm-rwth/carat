#include "typedef.h"
#include "tools.h"
#include "matrix.h"
#include "polyeder.h"

void vertex_standard(v)
vertex_TYP *v;
{
  int i,j;
  int w1;
 
  if(v->dim>0)
  {
    w1 = 0;
    i=0;
    while(i<v->dim && v->v[i] == 0)
       i++;
    if(i<v->dim)
      w1 = v->v[i];
    for(j=i+1;j<v->dim;j++)
    {
     if(v->v[j] != 0)
       w1 = GGT(w1, v->v[j]);
    }
    if(w1 < 0)
      w1 = -w1;
    if(w1 != 0)
    {
       for(j=i;j<v->dim;j++)
         v->v[j] /=w1;
    }
  }
}

void umsortieren(F,old_no, count, test, l)
fund_domain *F;
int old_no, *test, l;
int *count;
{
  int i;
  int d=0;
  for(i=0;i<old_no;i++)
  {
    if(count[i] < 0)
      free_vertex_fuber(&(F->vert[i]));
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
      free_wall_fuber(&(F->wall[i]));
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


void renumerate(v,test,l)
vertex_TYP *v;
int *test,l;
{
  int i,a;
  for(i=0;i<v->wall_no;i++)
  {
     a = v->wall[i];
     v->wall[i] = test[a];
  }
}


void wallAdd_vertex(i, v)
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


void streichen(v, test)
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
  

static int is_element(v, w)
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



void wallAdd_fund_domain(F, w)
fund_domain *F;
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
    F->wall = (wall_TYP **)realloc( (int *)F->wall,F->wall_SIZE *sizeof(wall_TYP *));
  }
  F->wall[F->wall_no] = w;
  F->wall_no++;
}



void vertexAdd_fund_domain(F, w)
fund_domain *F;
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
    F->vert = (vertex_TYP **)realloc( (int *)F->vert,F->vert_SIZE *sizeof(vertex_TYP *));
  }
  F->vert[F->vert_no] = w;
  F->vert_no++;
}


vertex_TYP *gis_neighbour(i,j,F,old_no)
int i,j;
fund_domain *F;
int old_no;
{
  int k,l,m;
  int a, u, tester1;
  int anz;
  vertex_TYP *erg;
  matrix_TYP *A;

  extern vertex_TYP *init_vertex_fuber();
  extern matrix_TYP *init_mat();
  extern int row_gauss();

  k=F->vert[j]->wall_no;
  if(k>F->vert[i]->wall_no)
    k=F->vert[i]->wall_no;
  erg = init_vertex_fuber(F->vert[0]->dim, k+1);
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
   { free_vertex_fuber(&erg); erg=NULL; return(NULL);}

  A = init_mat(anz,F->vert[0]->dim, "");
  for(k=0;k<anz;k++)
    for(l=0; l<F->vert[0]->dim;l++)
      A->array.SZ[k][l] = F->wall[erg->wall[k]]->gl[l];
  a = row_gauss(A);
  free_mat(A);
  if(a < F->vert[0]->dim-2)
   { free_vertex_fuber(&erg); erg=NULL; return(NULL);}
  return(erg);
}



vertex_TYP *is_neighbour(i,j,F,old_no)
int i,j;
fund_domain *F;
int old_no;
{
  int k,l,m;
  int a, u, tester1;
  int anz;
  vertex_TYP *erg;

  extern vertex_TYP *init_vertex_fuber();

  k=F->vert[j]->wall_no;
  if(k>F->vert[i]->wall_no)
    k=F->vert[i]->wall_no;
  erg = init_vertex_fuber(F->vert[0]->dim, k+1);
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
   { free_vertex_fuber(&erg); erg=NULL; return(NULL);}

  for(k=0; k<old_no ;k++)
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
       { free_vertex_fuber(&erg); erg=NULL; return(NULL);}
   }
  }
  return(erg);
}



int refine_fund_domain(F, h)
fund_domain *F;
wall_TYP *h;
{
 int i,j,k;
 int l;
 int old_no;
 int *count;
 int waste;
 int tester, *test;
 int anz;
 int p,n;
 vertex_TYP *v;

 extern void umsortieren();
 extern void renumerate();
 extern void streichen();
 extern int wall_times_vertex_fuber();
 extern int containes();
 extern void  wallAdd_fund_domain();
 extern void  vertexAdd_fund_domain();
 extern vertex_TYP *is_neighbour();
 

 p=0;n=0;
 waste = 0;
 old_no = F->vert_no;
 count = (int *)malloc(F->vert_no *sizeof(int *));
 for(i=0;i<F->vert_no;i++)
  count[i] = 0;
 for(i=0;i<old_no;i++)
 {
     count[i] = wall_times_vertex_fuber(h,F->vert[i]);
    if(count[i] > 0)
      p++;
    if(count[i] < 0)
      n++;
 }
 if(n == 0)
 {
   free(count);
   return(0);
 }
 for(i=0;i<old_no;i++)
 {
  if(count[i] < 0)
  {  
    for(j=0;j<old_no;j++)
    {
     if(count[j] > 0)
     {
        v = is_neighbour(i, j, F, old_no);
        if(v != NULL)
        {
           for(k=0; k<v->dim; k++)
           {
             v->v[k] = count[j] * F->vert[i]->v[k];
             v->v[k] -= count[i] * F->vert[j]->v[k];
           }
           vertex_standard(v);
           vertexAdd_fund_domain(F, v);
        }
     }
    }
  }
 }
  wallAdd_fund_domain(F, h);
  test = (int *)malloc((F->wall_no) *sizeof(int));
  l=0;
  for(i=0; i<F->wall_no-1;i++)
  {
    tester = 0;
    for(j=0;j<old_no && tester == 0; j++)
    {
     if(count[j] > 0)
       if(is_element(F->vert[j], i) == TRUE)
          tester = 1;
    }
    if(tester == 0)
     test[i] = -1;
    else
     {test[i] = l; l++;}
  }
  test[F->wall_no-1] = l;
  l++;

  for(i=0;i<old_no; i++)
  {
    if(count[i] == 0)
      wallAdd_vertex(F->wall_no-1, F->vert[i]);
  }

  for(i=0;i<old_no;i++)
  {
     if(count[i] ==0)
       streichen(F->vert[i],test );
     if(count[i] >=0)
       renumerate(F->vert[i],test,l);
  }
  for(i=old_no; i<F->vert_no;i++)
  {
     F->vert[i]->wall[F->vert[i]->wall_no] = F->wall_no-1;
     F->vert[i]->wall_no++;
     renumerate(F->vert[i],test,l);
  }

  umsortieren(F, old_no, count, test,l);
  
   free(count);
   free(test);
  return(1);
}
/*{{{}}}*/
