
/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: fub.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

#include "typedef.h"
#include "polyeder.h"
#include "matrix.h"
#include "symm.h"
#include "getput.h"
#include "bravais.h"
#include "sort.h"
#include "presentation.h"
#include "tools.h"
#include "datei.h"

static void kand_sort(ele_TYP **kand, int *rel_wand, int *k_anz);
static ele_TYP **create_new_cand( ele_TYP **old_kand, int k_anz,
	matrix_TYP **P_list, matrix_TYP **P_inv, polyeder_TYP *Pol, int anz);
static int isinpoly(matrix_TYP *Vec, polyeder_TYP *Pol);
static matrix_TYP *affine_vector_add(matrix_TYP *A, matrix_TYP *B, int n);
static void transform( matrix_TYP **elm, polyeder_TYP *Pol, matrix_TYP **trans,
	matrix_TYP **P_inv);
static matrix_TYP *mkpt(matrix_TYP *M, int *vert);
void normal_bary(int **B, int anz, int dim);
static int **bary_list(polyeder_TYP *Pol);
static int find(matrix_TYP *test, int **List, int anz);
static int *wall_bary(polyeder_TYP *Pol, int w);
static int satisfies_extra(polyeder_TYP *Pol, int w, matrix_TYP *test);
static int is_in_P(matrix_TYP *Vec, polyeder_TYP *Pol, int **zero, int *z);
static void add_vertex_to_polyeder(polyeder_TYP *F, vertex_TYP *w);
static void add_wall_to_polyeder(polyeder_TYP *F, wall_TYP *w);
static void add_wall_to_vertex(int i, vertex_TYP *v);
static int check_vertices( polyeder_TYP *Pol, int a, int b, int **w_vert_nr, 
	int *w_v_anz, matrix_TYP *inv);
static int vert_is_on_wall(matrix_TYP *v, wall_TYP *w, polyeder_TYP *Pol);
static int wall_vertex(int *w, vertex_TYP *v);
static void exchange(polyeder_TYP *Pol, int a, int b);
static vertex_TYP *is_neighbour(int i, int j, polyeder_TYP *F, int vertex_no);
static void calc_neighbors(polyeder_TYP *Pol);
static void extra_walls(polyeder_TYP *Pol, int a, int *w);
static void divide(polyeder_TYP *Pol, int w1, int w2, matrix_TYP *mat);
static void facts(polyeder_TYP *Pol, int **w_v_anz, int ***w_vert_nr);
static int check_polyeder(polyeder_TYP *Pol);
static void change_side_trafo(int i, polyeder_TYP *Pol, int k, 
	matrix_TYP **G_gen, matrix_TYP **G_inv, int G_gen_no, matrix_TYP *test);
static int look_for_wall(polyeder_TYP *Pol, int w, int **zero, int *anz, 
	int **w_vert_nr, int *w_v_anz, int *incase);
static void unterscheide(int i, int w1, polyeder_TYP *Pol, int *w_v_anz, 
	int **w_vert_nr, matrix_TYP *inv, int *incase);
static void check_side_trafo(polyeder_TYP *Pol, matrix_TYP **G_gen, int G_gen_no,
	 matrix_TYP **G_inv, matrix_TYP *Y);
static void Trans_anwenden( polyeder_TYP *Pol, int length, ele_TYP **orb,
	bravais_TYP *G, matrix_TYP **G_inv, matrix_TYP *Form, int P_anz,
	matrix_TYP *M, matrix_TYP **P_inv, matrix_TYP **P_list);


static void kand_sort(ele_TYP **kand, int *rel_wand, int *k_anz)
{
  int i;
  int s;

  s = 1;
  for(i=1; i< *k_anz; i++){
     if(rel_wand[i] == 0)
       free_ele(kand[i]);
     else{
     kand[s] = kand[i];
     s++;
     } 
  }
  kand = (ele_TYP **)realloc(kand, s*sizeof(ele_TYP *));
 *k_anz = s;
}/** kand_sort() **/

static ele_TYP **create_new_cand(old_kand,k_anz,P_list,P_inv,Pol,anz) 
ele_TYP **old_kand;
int k_anz;
matrix_TYP **P_list, **P_inv;
polyeder_TYP *Pol;
int anz;
{
  int		i,j,s;
  int		neu_anz;
  ele_TYP 	**kand;

 neu_anz = anz * (k_anz-1);
 kand = (ele_TYP **)malloc((neu_anz+1) * sizeof(ele_TYP*));
 for(i=0; i<neu_anz+1; i++){
    kand[i] = (ele_TYP *)malloc(sizeof(ele_TYP));
 }
 kand[0]->elm = copy_mat(old_kand[0]->elm);
 kand[0]->trans = NULL;
 kand[0]->schreier_vec = (int*)malloc(sizeof(int));
 kand[0]->schreier_vec[0] = 0;

 s = 1;
 for(i=1; i<k_anz; i++){
    for(j=0; j<anz; j++){
       kand[s]->elm = mat_mul(P_list[j],old_kand[i]->elm);
       Check_mat(kand[s]->elm);

       kand[s]->schreier_vec = 
		(int*)malloc((old_kand[i]->schreier_vec[0]+1) * sizeof(int));

       memcpy(kand[s]->schreier_vec, old_kand[i]->schreier_vec, 
				(old_kand[i]->schreier_vec[0]+1)*sizeof(int));

       kand[s]->trans = mat_mul(P_list[j],old_kand[i]->trans);
       Check_mat(kand[s]->trans);
/* 7.8.  kand[s]->trans = mat_mul(Pol->wall[j]->mat,old_kand[i]->trans); anne*/

       s++;
    }
 }
	if(s != neu_anz+1){
 		printf("fehler in create_neu_kand \n");
  		exit(3);
	}

return(kand);
}/** create_new_kand() **/

static int isinpoly(matrix_TYP *Vec, polyeder_TYP *Pol)
{
  int 	i, j, s;
  matrix_TYP 	*test,
		*M;
  
  j = Vec->rows;
  s = 0;

  M = init_mat(Pol->wall_no, j, "");

  for(i=0; i<M->rows; i++)
      memcpy(M->array.SZ[i], Pol->wall[i]->gl, j* sizeof(int));
  
  test = mat_mul(M,Vec);
  Check_mat(test);
 
  for (i=0; i<test->rows ; i++){
      if(test->array.SZ[i][0] < 0){
        free_mat(test); test = NULL;
        free_mat(M); M=NULL;
	return(i);			/* i = Nr der Wandgl. die verletzt ist*/
      }
  }
  free_mat(test); test = NULL;
  free_mat(M); M=NULL;
  return(-1); 				/*ist im Innern oder auf dem Rand*/ 
}/** isinpoly() **/

/* SPALTEN-KONVENTION, n=1 => A+B, n=-1 => A-B */
static matrix_TYP *affine_vector_add(A,B,n) 
matrix_TYP *A;
matrix_TYP *B;
int	    n;
{  
  matrix_TYP *C;
  int i,j;

  C = init_mat(A->rows, A->cols,"");
  for(i=0; i<C->rows-1; i++){
     for(j=0; j<C->cols; j++){
        if(n ==1)
        C->array.SZ[i][j] = 
		B->kgv*A->array.SZ[i][j] + A->kgv*B->array.SZ[i][j];
        else 
        C->array.SZ[i][j] = 
		B->kgv*A->array.SZ[i][j] - A->kgv*B->array.SZ[i][j];
     }
  }
  C->kgv = A->kgv * B->kgv;
  i = C->rows; 
  j = C->cols;
  C->array.SZ[i-1][j-1] = C->kgv;
  Check_mat(C);
return(C);
}/***affine_vector_add()***/

/*transformiere das Bahnelement elm in den Dirichletbereich D(x,Tx) durch 
  sukzessives Anwenden der Seitentransformationen, merke mir dies in trans */
static void transform(elm, Pol, trans, P_inv)
matrix_TYP **elm;
polyeder_TYP *Pol;
matrix_TYP **trans;
matrix_TYP **P_inv;
{
  int		test;
  matrix_TYP	*tmp,*tmp1, *tmp2;
  matrix_TYP	*erg;

  erg = copy_mat(elm[0]);
  test = isinpoly(erg, Pol);
  while(test != -1){
     tmp = mat_mul(P_inv[test], erg);
     Check_mat(tmp);
     free_mat(erg); erg = NULL;
     erg = tmp;
     test = isinpoly(erg, Pol);
  }
  tmp1 = affine_vector_add(erg,elm[0],-1);  /* gx+t = g'x = erg */
  
  if(trans[0] != NULL){
  tmp2 = affine_vector_add(trans[0],tmp1,1); /*addiere tmp1 zu trans hinzu*/
  free_mat(trans[0]); trans[0] = NULL;
  trans[0] = tmp2; 
  free_mat(tmp1), tmp1 = NULL;
  }
  else
  trans[0] = tmp1; 
  free_mat(elm[0]); elm[0] = NULL;
  elm[0] = erg;
}/** transform() **/

static matrix_TYP *mkpt(matrix_TYP *M, int *vert)
{
  int		i;
  int 		vkgv, kkg;
  matrix_TYP 	*Y;

  vkgv = vert[M->rows-1];
  Y = init_mat(M->rows,M->cols,"");
  if(M->kgv == vert[M->rows-1]){
     for(i=0; i<Y->rows-1; i++)
        Y->array.SZ[i][0] = (M->array.SZ[i][0]) + vert[i];
     Y->array.SZ[Y->rows-1][0] = 2*(M->kgv);  
     Y->kgv = 2*(M->kgv);
     Check_mat(Y);
  }
  else{
     kkg = 2*(M->kgv)*vkgv;
     for(i=0; i<Y->rows-1; i++)
        Y->array.SZ[i][0] = 
		vkgv*(M->array.SZ[i][0]) + (M->kgv)*(vert[i]);
     Y->array.SZ[Y->rows-1][0] = kkg;
     Y->kgv = kkg;
     Check_mat(Y);
  }
return(Y);
}/** mkpt() **/

void normal_bary(B,anz,dim)
int **B;
int anz;
int dim;
{
  int i,j,l;
  int w1,w2;

  for(l=0; l<anz; l++){ 
     if(dim>0)
     {
       i=0;
       while(i<dim && B[l][i] == 0)
          i++;
       if(i<dim)
         w1 = B[l][i];
       for(j=i+1;j<dim && w1 != 1 && w1 != -1;j++)
       {
        if(B[l][j] != 0)
        {
          w2 = GGT(w1, B[l][j]);
          w1 = w2;
        }
       }
       if(w1 < 0)
         w1 = -w1;
       if(w1 != 0)
       {
          for(j=0;j<dim;j++)
            B[l][j] /= w1;
       }
     }
  }/* for(l) */
}/* normal_bary() */

/*berechne baryzentrum der Waende von P */
static int **bary_list(polyeder_TYP *Pol)
{
  int **erg;
  int *v_wand;	 	/*Nummern der Ecken der Wand*/
  int v_anz; 		/*Anzahl der Ecken der Wand*/
  int i, j, k, a, b;
  int *list;
  int kgv;
  int d;

  d = Pol->vert[0]->dim;
  erg = (int **)malloc(Pol->wall_no * sizeof(int*));
  for(i=0; i<Pol->wall_no; i++){
     erg[i] = (int *)malloc(d * sizeof(int));
     for(k=0; k<d; k++)
         erg[i][k] = 0;
     }

  for(k=0; k<Pol->wall_no; k++){

  	/* bestimme v_wand und v_anz */
        v_anz = 0;
	v_wand = (int *) calloc(Pol->vert_no , sizeof(int));  

	/* bestimme Eckenmenge der Wand*/
	for(a=0; a<Pol->vert_no; a++)     
	{ 
	   for(b=0; b<Pol->vert[a]->wall_no; b++)
	    {
	     if(Pol->vert[a]->wall[b] == k)
             {	
	      	v_wand[v_anz] = a;
              	v_anz++;
             }
            }
 	}
	v_wand = (int*)realloc(v_wand,v_anz * sizeof(int));

	/*bestimme Baryzentrum der Wand*/
      list = (int*) calloc(v_anz ,sizeof(int));
      for(i=0; i<v_anz; i++) 		
      {
	/*v_wand[i]= Originalnummer der i-tn Ecke im Ausgangspolyeder*/
      list[i] = Pol->vert[v_wand[i]]->v[d-1];
      }
      kgv = KKGV(list, v_anz);	/*funktion KKGV extern */ 
  
      for(j=0; j<d; j++){
          for(i=0; i<v_anz; i++)
          {
          	erg[k][j] +=  (kgv/list[i])*(Pol->vert[v_wand[i]]->v[j]); 
    	    /* teile noch durch letzte Komponente, diese ist Anzahl*kgv */
          }
      }
    free(list);
   free(v_wand);
   }
normal_bary(erg,Pol->wall_no,Pol->vert[0]->dim);
return(erg);
}/**bary_list()**/

/*test ist Spalte*/
static int find(matrix_TYP *test, int **List, int anz)
{
  int	i,k,s,l;
  
  k = 0;
  for(i=0; i<anz; i++){
      for(s=0; s<test->rows; s++){
      	  if(test->array.SZ[s][0] == List[i][s])
          k++;
      }
      l=(k == test->rows);
  if(l != 0)
  return(i);
  k = 0;
  }
 return(-1);
}/**find()**/

static int *wall_bary(polyeder_TYP *Pol, int w)
{
  int *erg;
  int *v_wand;	 	/*Nummern der Ecken der Wand*/
  int v_anz; 		/*Anzahl der Ecken der Wand*/
  int i, j, a, b;
  int *list;
  int kgv;
  int d;

  d = Pol->vert[0]->dim;
  erg = (int *)calloc(d , sizeof(int));

  	/* bestimme v_wand und v_anz */
        v_anz = 0;
	v_wand = (int *) calloc(Pol->vert_no , sizeof(int));  

	/* bestimme Eckenmenge der Wand*/
	for(a=0; a<Pol->vert_no; a++)     
	{ 
	   for(b=0; b<Pol->vert[a]->wall_no; b++)
	    {
	     if(Pol->vert[a]->wall[b] == w)
             {	
	      	v_wand[v_anz] = a;
              	v_anz++;
             }
            }
 	 }
	 v_wand = (int*)realloc(v_wand,v_anz * sizeof(int));

	/*bestimme Baryzentrum der Wand*/
      list = (int*) calloc(v_anz ,sizeof(int));
      for(i=0; i<v_anz; i++) 		
      {
	/*v_wand[i]= Originalnummer der i-tn Ecke im Ausgangspolyeder*/
      list[i] = Pol->vert[v_wand[i]]->v[d-1];
      }
      kgv = KKGV(list, v_anz);	/*funktion KKGV extern */ 
  
      for(j=0; j<d; j++){
          for(i=0; i<v_anz; i++)
          {
          	erg[j] +=  (kgv/list[i])*(Pol->vert[v_wand[i]]->v[j]); 
    	    /* teile noch durch letzte Komponente, diese ist Anzahl*kgv */
          }
      }
    free(list);
   
   free(v_wand);

return(erg);
}/** wall_bary() **/

static int satisfies_extra(polyeder_TYP *Pol, int w, matrix_TYP *test)
{ 
  int i,j,e;
  
  for(i=0; i<Pol->wall[w]->ext_no; i++){
     e = 0;
     for(j=0; j<Pol->vert[0]->dim; j++)
        e += Pol->wall[w]->extra[i][j] * test->array.SZ[j][0];
     if(e < 0)
       return(-1);
  }
  for(i=0; i<Pol->wall[w]->next_no; i++){
     e = 0;
     for(j=0; j<Pol->vert[0]->dim; j++)
        e += Pol->wall[w]->next[i][j] * test->array.SZ[j][0];
     if(e < 0)
       return(-1);
  }
return(1);
}/** satisfies_extra() **/

static int is_in_P(matrix_TYP *Vec, polyeder_TYP *Pol, int **zero, int *z)
{
  int 	i, j, s;
  int   z1;
  matrix_TYP 	*test,
		*M;
  
  j = Vec->rows;
  s = 0;
  z1 = 0;

  zero[0] = (int*)calloc(Pol->wall_no, sizeof(int));
  M = init_mat(Pol->wall_no, j, "");

  for(i=0; i<M->rows; i++){
      memcpy(M->array.SZ[i], Pol->wall[i]->gl, j* sizeof(int));
  }

  test = mat_mul(M,Vec);
  Check_mat(test);
 
  for (i=0; i<Pol->wall_no ; i++){
      if(test->array.SZ[i][0] < 0){
        free_mat(M); M = NULL;
        free_mat(test); test = NULL;
 	s = (-1)* (i+1);
	return(s);
      }
  }
  for (i=0; i<Pol->wall_no ; i++){
      if(test->array.SZ[i][0] == 0){
           s = satisfies_extra(Pol,i,Vec);
           if(s == 1) {
              zero[0][z1++] = i;
           }
      }
  }
  *z = z1;
  if(z1 == 1){
    free_mat(M); M = NULL;
    free_mat(test); test = NULL;
    return(zero[0][0]+1); 
  }
  free_mat(M); M = NULL;
  free_mat(test); test = NULL;
  return(0); 
}/** is_in_P() **/

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

static int check_vertices(Pol,a,b,w_vert_nr, w_v_anz,inv)
polyeder_TYP *Pol;
int a,b;		/* wand a wird auf wand b abgeb. via (w(a)->mat)^(-1) */
int **w_vert_nr;
int *w_v_anz;
matrix_TYP *inv;
{
 int	i,j,k,s,t;
 int	test;
 matrix_TYP *tmp1, *tmp2;

 tmp1 = init_mat(Pol->vert[0]->dim, 1,"");
 for(i=0; i<w_v_anz[a]; i++){
    k = w_vert_nr[a][i];
    for(j=0; j<Pol->vert[0]->dim; j++)
       tmp1->array.SZ[j][0] = Pol->vert[k]->v[j];
    tmp2 = mat_mul(inv,tmp1);
    Check_mat(tmp2);
    test = 0;
    for(j=0; j<w_v_anz[b]; j++){
       k = w_vert_nr[b][j];
       s = 0;
       for(t=0; t<Pol->vert[0]->dim; t++) {
          if( tmp2->array.SZ[t][0] == Pol->vert[k]->v[t])
               s++;
       }
       if(s == Pol->vert[0]->dim){
         test = 1;
         break;				/* 21.11.97 */
       }
    }  
    if(test == 0){
      free_mat(tmp1); tmp1 = NULL;
      free_mat(tmp2); tmp2 = NULL;
      return(i);	/*Nr. der ersten Ecke die nicht gut abgebildet wird*/
    }
    free_mat(tmp2); tmp2 = NULL;
 }/** for(i) **/
 free_mat(tmp1); tmp1 = NULL;
 if(w_v_anz[a] == w_v_anz[b])
   return(-1);
 else{
   printf(" w_v_anz[a] != w_v_anz[b] \n");
   return(-3);
 }
}/** check_vertices() **/

static int vert_is_on_wall(v,w,Pol)
matrix_TYP *v;
wall_TYP *w;
polyeder_TYP *Pol;
{
  int i,j,s;
 
  j = 0;
  for(i=0; i<Pol->vert[0]->dim; i++)
      j += w->gl[i] * (v->array.SZ[i][0]);
  if(j != 0) 
    return(-1);
  else{
    for(s=0; s<w->next_no; s++){
        j = 0;
        for(i=0; i<Pol->vert[0]->dim; i++)
            j += w->next[s][i] * (v->array.SZ[i][0]);
        if( j < 0 )
          return(-1);
    }
    for(s=0; s<w->ext_no; s++){
        j = 0;
        for(i=0; i<Pol->vert[0]->dim; i++)
            j += w->extra[s][i] * (v->array.SZ[i][0]);
        if( j < 0 )
          return(-1);
    }
  } 
return(1);
}/** vert_is_on_wall() **/

static int wall_vertex(int *w, vertex_TYP *v)
{
  int i,e;
  
  e = 0;
  for(i=0; i<v->dim; i++)
      e += w[i] * v->v[i];

  return(e);
}/** wall_vertex() **/

static void exchange(polyeder_TYP *Pol, int a, int b)
{
  int 	i,s;
 
  s = 0;
  for(i=0; i<Pol->vert[b]->wall_no; i++){
     if(a == Pol->vert[b]->wall[i]){
        s = i;
        i = Pol->vert[b]->wall_no-1;
     }
  }
  memmove(Pol->vert[b]->wall+s, Pol->vert[b]->wall+s+1,
                                      (Pol->vert[b]->wall_no - s)*sizeof(int));
  Pol->vert[b]->wall[i-1] = Pol->wall_no-1;

}/** exchange() **/

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

static void calc_neighbors(Pol)
polyeder_TYP *Pol;
{
  int		i,j,k,l,s;
  int		dim;
  int		v_anz;
  int		w_anz;
  int		*count;
  int		*w_v_anz;
  int		**w_vert_nr;
  int		**v_wand_nr;
 
  dim = Pol->vert[0]->dim;
  v_anz = Pol->vert_no;
  w_anz = Pol->wall_no;

  count = (int*)calloc(w_anz, sizeof(int));
  w_v_anz = (int*)malloc(w_anz*sizeof(int));
  w_vert_nr = (int**)malloc(w_anz * sizeof(int*));
  v_wand_nr = (int**)malloc(v_anz * sizeof(int*));

  for(i=0; i<w_anz; i++)
  {
    s = 0;
    w_vert_nr[i] = (int*)malloc(v_anz *sizeof(int));
    for(j=0; j<v_anz; j++)
       for(k=0; k<Pol->vert[j]->wall_no; k++)
          if(Pol->vert[j]->wall[k] == i)
            w_vert_nr[i][s++] = j;
    w_v_anz[i] = s;
    w_vert_nr[i] = (int*)realloc(w_vert_nr[i], s * sizeof(int));
  }

  for(i=0; i<v_anz; i++)
  {
     v_wand_nr[i] = (int*)calloc(w_anz , sizeof(int));
     for(j=0; j<Pol->vert[i]->wall_no; j++)
        v_wand_nr[i][j] = Pol->vert[i]->wall[j];
  }

  for(i=0; i<w_anz; i++)
  {
    if(Pol->wall[i]->next == NULL)
      Pol->wall[i]->next = (int**)malloc(w_anz*sizeof(int*));
    else
      Pol->wall[i]->next = (int**)realloc(Pol->wall[i]->next, w_anz * sizeof(int*));
    s = 0;
    for(l=0; l<w_anz; l++){
       if(l!=i){
          count[l] = 0;
          for(k=0; k<w_v_anz[i]; k++){
             for(j=0; j<w_v_anz[l]; j++)
                if(w_vert_nr[i][k] == w_vert_nr[l][j])     
                   count[l] ++;
          }
          if(count[l] >= dim-2){
             Pol->wall[i]->next[s] = (int*)malloc(dim * sizeof(int));
             memcpy(Pol->wall[i]->next[s],Pol->wall[l]->gl, dim*sizeof(int));
             s++;
          }
       }
    }
    Pol->wall[i]->next = (int **)realloc(Pol->wall[i]->next,s*sizeof(int*)); 
    Pol->wall[i]->next_no = s; 
  }

free(count);
for(i=0; i<v_anz; i++)
    free(v_wand_nr[i]);
free(v_wand_nr);
for(i=0; i<w_anz; i++)
    free(w_vert_nr[i]);
free(w_vert_nr);
free(w_v_anz);
}/** calc_neighbors() **/

static void extra_walls(Pol,a,w)
polyeder_TYP *Pol;
int a;
int *w;
{
  int i;
  int b;
  int a_anz;
  int dim;

  dim  = Pol->vert[0]->dim;
  b = Pol->wall_no-1;
  a_anz = Pol->wall[a]->ext_no;

  if(Pol->wall[b]->extra == NULL)
     Pol->wall[b]->extra = (int **)malloc((a_anz+1) * sizeof(int*));
  else
     Pol->wall[b]->extra = (int **)realloc(Pol->wall[b]->extra,(a_anz+1) * sizeof(int*));
  if(Pol->wall[a]->extra == NULL)
     Pol->wall[a]->extra = (int **)malloc((a_anz+1) * sizeof(int*));
  else 
     Pol->wall[a]->extra = (int **)realloc(Pol->wall[a]->extra,(a_anz+1) * sizeof(int*));

     Pol->wall[b]->extra[a_anz]= (int*)malloc(dim * sizeof(int));
     Pol->wall[a]->extra[a_anz]= (int*)malloc(dim * sizeof(int));

  Pol->wall[a]->ext_no ++;
  for(i=0; i<dim; i++)
     Pol->wall[a]->extra[a_anz][i] = w[i];
  Pol->wall[b]->ext_no ++;
  for(i=0; i<dim; i++)
     Pol->wall[b]->extra[a_anz][i] = (-1)*w[i];
}/** extra_walls() **/

/* unterteile w1 mit dem Bild von w2 */
static void divide(polyeder_TYP *Pol, int w1, int w2, matrix_TYP *mat)
{
  int		i,j,l,ii,jj;
  int		s,k;
  int		a,b;
  int 		n,p,z;
  int 		nn,pp,zz;
  int		dim;
  int		next_no,ext_no,no;
  int		old_no;
  int		wall_length;
  int		wall_neu;
  int		a_v_anz;
  int		*w;
  int		*wall_list;
  int		*a_vert_nr;
  int		*count;
  int		*ccount;
  int 		*wall_test, *necessary;
  int		*v_on_a;
  int 		SIZE;
  vertex_TYP 	*v;
  matrix_TYP	*extra_gl;
  matrix_TYP	*new_gl;

  SIZE = 100;
  dim = Pol->vert[0]->dim;
  next_no = Pol->wall[w1]->next_no;
  ext_no = Pol->wall[w1]->ext_no;
  no = Pol->wall[w1]->ext_no+Pol->wall[w1]->next_no;
  extra_gl = init_mat(no, dim, "");
  wall_list = (int*)malloc(SIZE * sizeof(int));
  v = NULL;

  wall_list[0] = w2;
  wall_length = 1;

  for(i=0; i<ext_no; i++){
     memcpy(extra_gl->array.SZ[i],Pol->wall[w1]->extra[i],dim * sizeof(int));
  }
  for(i=ext_no; i<next_no; i++){
     memcpy(extra_gl->array.SZ[i],Pol->wall[w1]->next[i],dim * sizeof(int));
  }

  new_gl = mat_mul(extra_gl,mat);
  	/*   new_gl = mat_mul(extra_gl,Pol->wall[w1]->mat);   */
  Check_mat(new_gl);

  free_mat(extra_gl); extra_gl = NULL;

  w = (int*)calloc(dim, sizeof(int));

  for(i=0; i<no; i++){
      memcpy(w,new_gl->array.SZ[i], dim*sizeof(int));
      wall_neu = 0;

      for(j=0; j<wall_length; j++){
          a = wall_list[j]; 
          a_vert_nr = (int*)calloc(Pol->vert_no , sizeof(int));
          v_on_a = (int*)calloc(Pol->vert_no, sizeof(int));
          s = 0;
          for(jj=0; jj<Pol->vert_no; jj++){
             for(k=0; k<Pol->vert[jj]->wall_no; k++)
                if(Pol->vert[jj]->wall[k] == a){
                  a_vert_nr[s++] = jj;
                  v_on_a[jj] = 1;
                }
          }
          a_v_anz = s;

          a_vert_nr = (int *)realloc(a_vert_nr, s*sizeof(int));

          count = (int*)calloc(a_v_anz, sizeof(int));
          p = 0; n = 0; z = 0;
          pp = 0; nn = 0; zz = 0;

          for(l=0; l<a_v_anz; l++){
             count[l] = wall_vertex(w,Pol->vert[a_vert_nr[l]]);
             if(count[l] > 0)
                p++;
             if(count[l] < 0)
                n++;
             if(count[l] == 0)
                z++;
          } /** for(l<a_v_anz ) **/

          ccount = (int*)calloc(Pol->vert_no, sizeof(int));
          pp = 0; nn = 0; zz = 0;
          for(l=0; l<Pol->vert_no; l++){
             ccount[l] = wall_vertex(w,Pol->vert[l]);
             if(ccount[l] > 0)
                pp++;
             if(ccount[l] < 0)
                nn++;
             if(ccount[l] == 0)
                zz++;
          } /** for(l<Pol->vert_no) **/

          if(n!=0 && p!=0){
/* bestimme neue Eckenmenge */

             old_no = Pol->vert_no;
             wall_test = (int *)malloc(Pol->wall_no *sizeof(int));
             for(ii=0;ii<Pol->wall_no;ii++)
                 wall_test[ii] = 2;
             for(ii=0;ii<old_no;ii++)
             {
               if(ccount[ii]>0)
               {
                  for(jj=0;jj<Pol->vert[ii]->wall_no;jj++)
                  {
                    if(wall_test[Pol->vert[ii]->wall[jj]] == 2)
                      wall_test[Pol->vert[ii]->wall[jj]] = 1;
                    if(wall_test[Pol->vert[ii]->wall[jj]] == -1)
                      wall_test[Pol->vert[ii]->wall[jj]] = 0;
                  }
               }
               if(ccount[ii] < 0)
               {
                  for(jj=0;jj<Pol->vert[ii]->wall_no;jj++)
                  {
                    if(wall_test[Pol->vert[ii]->wall[jj]] == 2)
                      wall_test[Pol->vert[ii]->wall[jj]] = -1;
                    if(wall_test[Pol->vert[ii]->wall[jj]] == 1)
                      wall_test[Pol->vert[ii]->wall[jj]] = 0;
                  }
               }
             }
             necessary = (int *)malloc(Pol->vert_no *sizeof(int));
             for(ii=0;ii<old_no;ii++)
               necessary[ii] = FALSE;
             for(ii=0;ii<old_no;ii++)
             {
              if(ccount[ii] != 0)
              {
               k=0;
               for(jj=0;jj<Pol->vert[ii]->wall_no && necessary[ii] == FALSE;jj++)
               {
                 if( wall_test[Pol->vert[ii]->wall[jj]] == 0)
                   k++;
                 if(k>= Pol->vert[ii]->dim-2)
                 {
                   if(ccount[ii] > 0)
                     necessary[ii] = 1;
                   if(ccount[ii] < 0)
                     necessary[ii] = -1;
                   if(ccount[ii] == 0)
                     necessary[ii] = 0;
                 }
               }
              }
             }
            

             if(Pol->vert[0]->dim == 2)
             {
               for(ii=0;ii<old_no;ii++)
               {
                 if(ccount[ii] > 0)
                   necessary[ii] = 1;
                 if(ccount[ii] < 0)
                   necessary[ii] = -1;
                 if(ccount[ii] == 0)
                   necessary[ii] = 0;
               }
             }
            
            
             for(ii=0;ii<old_no;ii++)
             {
              if(necessary[ii] == -1)
              {  
                for(jj=0;jj<old_no;jj++)
                {
                 if(necessary[jj] == 1)
                 {
                    if(v_on_a[ii]==1 && v_on_a[jj]==1)
                      v = is_neighbour(ii, jj, Pol, old_no);
                    if(v != NULL)
                    {
                      for(k=0; k<v->dim; k++)
                        v->v[k] = ccount[jj]* Pol->vert[ii]->v[k] 
                                  - ccount[ii]* Pol->vert[jj]->v[k];
                      normal_vertex(v);
                      add_vertex_to_polyeder(Pol, v);
                    }
                    v = NULL;
                 }
                }
              }
             }

            for(ii=old_no; ii < Pol->vert_no ;ii++)
            {
               if(Pol->vert[ii]->wall_no == Pol->vert[ii]->wall_SIZE)
               {
                 Pol->vert[ii]->wall_SIZE += extsize1;
                 Pol->vert[ii]->wall = (int *)realloc(Pol->vert[ii]->wall,
                                      Pol->vert[ii]->wall_SIZE *sizeof(int));
               }
               Pol->vert[ii]->wall[Pol->vert[ii]->wall_no] = Pol->wall_no;
               Pol->vert[ii]->wall_no ++; 
            }/** for(old_no<ii<Pol->vert_no) **/

             free(necessary);
             free(wall_test);

             add_wall_to_polyeder(Pol,Pol->wall[a]);
             extra_walls(Pol,a,w);

             for(l=0; l<a_v_anz; l++){
                 b = a_vert_nr[l];
                 if(count[l] == 0)
                    add_wall_to_vertex(Pol->wall_no-1, Pol->vert[b]);
                 if(count[l]<0)
                    exchange(Pol,a,b);
             } /** for(l<a_v_anz) **/

             if(wall_length+wall_neu == SIZE){
                SIZE += SIZE;
                wall_list = (int *)realloc(wall_list, SIZE * sizeof(int)); 
             }
             wall_list[wall_length+wall_neu] = Pol->wall_no-1;
             wall_neu ++;
      
          } /** if(n!=0&&p!=0) **/

          free(count);
          free(ccount);
          free(a_vert_nr);
          free(v_on_a);
      } /** for(j<wall_length) **/
             
         wall_length += wall_neu;

  } /** for(i<no = next_no+ext_no) **/

  free(wall_list);
  free(w);

} /* divide() */

static void facts(Pol,w_v_anz,w_vert_nr)
polyeder_TYP *Pol;
int **w_v_anz;
int ***w_vert_nr;
{
  int	l,s,j,k;

    w_v_anz[0] = (int*)malloc(Pol->wall_no *sizeof(int));
    w_vert_nr[0] = (int**)malloc(Pol->wall_no * sizeof(int*));

    for(l=0; l<Pol->wall_no; l++){
      s = 0;
      w_vert_nr[0][l] = (int*)malloc(Pol->vert_no *sizeof(int));
      for(j=0; j<Pol->vert_no; j++)
         for(k=0; k<Pol->vert[j]->wall_no; k++)
            if(Pol->vert[j]->wall[k] == l)
              w_vert_nr[0][l][s++] = j;
      w_v_anz[0][l] = s;
      w_vert_nr[0][l] = (int*)realloc(w_vert_nr[0][l], s * sizeof(int));
    }
}/** facts() **/

static int check_polyeder(polyeder_TYP *Pol)
{
  int i,j,k,ok;
  int **Bary;
  matrix_TYP *test,*inv,*tmp;

  test = init_mat(Pol->vert[0]->dim,1,"");
  Bary = bary_list(Pol);
  
  for(i=0; i<Pol->wall_no; i++){
    for(j=0; j<test->rows; j++)
       test->array.SZ[j][0] = Bary[i][j];
    inv = mat_inv(Pol->wall[i]->mat);
    tmp = mat_mul(inv,test);	
    Check_mat(tmp);
		/*tmp liegt auf dem Rand des vorigen Polyeders */
    ok = find(tmp,Bary, Pol->wall_no);
      if(ok == -1){
        free_mat(test); test = NULL;
        free_mat(tmp); tmp = NULL;
        free_mat(inv); inv = NULL;
        for(k=0; k<Pol->wall_no; k++)
           free(Bary[k]);
        free(Bary);
        return(i+1);
      }
    free_mat(tmp); tmp = NULL;
    free_mat(inv); inv = NULL;
  }
  free_mat(test); test = NULL;
  for(k=0; k<Pol->wall_no; k++)
     free(Bary[k]);
  free(Bary);
  return(0);
}/** check_polyeder () **/


static void change_side_trafo(i,Pol,k,G_gen,G_inv,G_gen_no,test)
int i;
polyeder_TYP *Pol;
int k;
matrix_TYP **G_gen;
matrix_TYP **G_inv;
int G_gen_no;
matrix_TYP *test;
{
 word_TYP 	*word;
 matrix_TYP 	*tmp, *tmp2;
 int		j,s;

   if(Pol->wall[k]->neu == 1){
      printf(" habe alte Wand ueberschritten \n");
      exit(3);
   }
 /* do{*/
   tmp2 = mat_mul(Pol->wall[i]->mat, Pol->wall[k]->mat); 
   Check_mat(tmp2);
   free_mat(Pol->wall[i]->mat); Pol->wall[i]->mat = NULL;
   Pol->wall[i]->mat = tmp2;
   word = init_word(Pol->vert[0]->dim);
   Word_mul(Pol->wall[i]->word, Pol->wall[k]->word,G_gen, G_inv, G_gen_no,
            word);
   free_word(Pol->wall[i]->word);
   Pol->wall[i]->word = word;
   tmp = mat_mul(Pol->wall[i]->mat,test);
   Check_mat(tmp);
   s = 0;
   for(j=0; j<Pol->vert[0]->dim; j++)
      s += Pol->wall[k]->gl[j] * tmp->array.SZ[j][0];

   free_mat(tmp); tmp = NULL;

/* }while(s < 0); */

}/** change_side_trafo() **/

static int look_for_wall(Pol,w,zero,anz,w_vert_nr, w_v_anz, incase)
polyeder_TYP *Pol;
int w;
int **zero;
int *anz;
int **w_vert_nr;
int *w_v_anz;
int *incase;
{
  int 		i,j;
  int		a,b,z;
  int		count;
  matrix_TYP 	*v, *tmp,*inv;

  v = init_mat(Pol->vert[0]->dim,1,"");

  inv = mat_inv(Pol->wall[w]->mat); Check_mat(inv);

  count = 0;
  for(i=0; i<anz[0]; i++){
     b = 0;
     z = zero[0][i];
     for(j=0; j<w_v_anz[w]; j++){
           for(a=0; a<v->rows; a++)
              v->array.SZ[a][0] = Pol->vert[w_vert_nr[w][j]]->v[a];
           tmp = mat_mul(inv,v); Check_mat(tmp);
           for(a=0; a<v->rows; a++)
               b += Pol->wall[z]->gl[a] * tmp->array.SZ[a][0]; 
        if(b != 0)
          break;
     } 
     if(b == 0){
       zero[0][count ++] = z;
     }
  }
  
  if(count == 1){ 
       anz[0] = 1;
       free_mat(v); v = NULL;
       *incase = 1;
       return(zero[0][0]);
  }
/*  else  */
  if(count > 1){
       anz[0] = count;
       zero[0] = (int*)realloc(zero[0], count * sizeof(int));

       b = 0;
       for(i=0; i<anz[0]; i++){
          z = zero[0][i];
          count = 0;
          for(j=0; j<w_v_anz[w]; j++){
              for(a=0; a<v->rows; a++)
                 v->array.SZ[a][0] = Pol->vert[w_vert_nr[w][j]]->v[a];
              tmp = mat_mul(inv,v); Check_mat(tmp);
              a = vert_is_on_wall(tmp,Pol->wall[z],Pol);
              if(a == -1)
                 break;
              count ++;
              if(count == w_v_anz[w]){
                 b++;
              }
          } 
       }
       if(b == 1){
              *incase = 2;
              free_mat(v); v = NULL;
              return(zero[0][i]);
       }
       *incase = 3;
       free_mat(v); v = NULL;
       return(-1);
  }

 printf(" war wohl nichts \n");
 exit(3); 
}/** look_for_wall() **/

static void unterscheide(i,w1,Pol,w_v_anz,w_vert_nr,inv,incase)
int i,w1;
polyeder_TYP *Pol;
int *w_v_anz;
int **w_vert_nr;
matrix_TYP *inv;
int *incase;
{
  int 		s,s1;
  int		j,k;
  matrix_TYP	*tmp1, *tmp2;

          tmp1 = init_mat(Pol->vert[0]->dim,1,"");

          s = 0;
	  s1 = 0;
	  for(j=0; j<w_v_anz[w1]; j++){
	     for(k=0; k<Pol->vert[0]->dim; k++)
	        tmp1->array.SZ[k][0] = Pol->vert[w_vert_nr[w1][j]]->v[k];
     /*	
     tmp2 = mat_mul(inv,tmp1);
     */    
	     tmp2 = mat_mul(Pol->wall[i]->mat,tmp1);
	     Check_mat(tmp2);
	     s1 = vert_is_on_wall(tmp2,Pol->wall[i],Pol);
	     free_mat(tmp2); tmp2 = NULL;
	     if(s1 == -1)
	     break;
	  }

          for(j=0; j<w_v_anz[i]; j++){
             for(k=0; k<Pol->vert[0]->dim; k++)
                tmp1->array.SZ[k][0] = Pol->vert[w_vert_nr[i][j]]->v[k];
             tmp2 = mat_mul(inv,tmp1);
             Check_mat(tmp2);
             s = vert_is_on_wall(tmp2,Pol->wall[w1],Pol);
             free_mat(tmp2); tmp2 = NULL;
             if(s == -1)
             break;
          }
          free_mat(tmp1); tmp1 = NULL;
if(s != -1 && s1 == -1)
   *incase = 2;
if((s == -1 && s1 != -1) || (s == -1 && s1 == -1))
   *incase = 1;
}/** unterscheide () **/

static void check_side_trafo(Pol,G_gen, G_gen_no,G_inv, Y)
polyeder_TYP *Pol;
matrix_TYP **G_gen;
int G_gen_no;
matrix_TYP **G_inv;
matrix_TYP *Y;
{
 int 		i,j,k,l;
 int		a;
 int		s,t;
 int		z;
 int		isok;
 int		old_P_wall_no;
 int 		incase;
 int		dim;
 int		*inv_w;
 int		**Bary;
 int		*bar;
 int		*check;
 int 		*zero;
 int		*w_v_anz;
 int		**w_vert_nr;
 matrix_TYP 	*tmp, *A, *B;
 matrix_TYP 	*test;
 matrix_TYP 	*inv;
 
 extern void Word_mul();
 
 incase = 0;
 dim = Pol->vert[0]->dim;

 Bary = bary_list(Pol);
 test = init_mat(Pol->vert[0]->dim,1,"");
 check = (int*)calloc(Pol->wall_no, sizeof(int));

 calc_neighbors(Pol);

 w_v_anz = (int*)malloc(Pol->wall_no *sizeof(int));
 w_vert_nr = (int**)malloc(Pol->wall_no * sizeof(int*));

 for(l=0; l<Pol->wall_no; l++){
   s = 0;
   w_vert_nr[l] = (int*)malloc(Pol->vert_no *sizeof(int));
   for(j=0; j<Pol->vert_no; j++)
     for(k=0; k<Pol->vert[j]->wall_no; k++)
       if(Pol->vert[j]->wall[k] == l)
	 w_vert_nr[l][s++] = j;
   w_v_anz[l] = s;
   w_vert_nr[l] = (int*)realloc(w_vert_nr[l], s * sizeof(int));
 }

 i = 0;
 do{
   do{
     bar = wall_bary(Pol,i);
     for(j=0; j<test->rows; j++)
       test->array.SZ[j][0] = bar[j];
     inv = mat_inv(Pol->wall[i]->mat);
     tmp = mat_mul(inv,test);	
     Check_mat(tmp);
     /*tmp liegt auf dem Rand des vorigen Polyeders */

     t = is_in_P(tmp,Pol,&zero,&z);
     if(t < 0) /* tmp nicht in Pol, wand Nr abs(t)-1 ist verletzt. */
       { 
         /*aendere die Seitentrafo zu i-ter wand und gehe damit nochmal in die             Schleife*/
	 k = abs(t)-1;
	 change_side_trafo(i,Pol,k,G_gen,G_inv,G_gen_no,test);
	 
       }
     if(t>0){  /* tmp auf Rand von Pol  und zwar auf Wand t-1*/
       isok = check_vertices(Pol,i,t-1,w_vert_nr, w_v_anz,inv);
       if(isok == -1){
	 check[i] = 1; 
	 if(i != t-1){
	   check[t-1] = 1;
	   if(mat_comp(inv,Pol->wall[t-1]->mat) != 0){
	     free_mat(Pol->wall[t-1]->mat);
	     Pol->wall[t-1]->mat = copy_mat(inv);
	     if(Pol->wall[i]->word->word != NULL){
	       /* invertiere Wort !!! */
	       inv_w = (int*)calloc(Pol->wall[i]->word->word[0] +1, sizeof(int));
	       inv_w[0] = Pol->wall[i]->word->word[0];
	       for(a = 1; a <= Pol->wall[i]->word->word[0]; a++)
		 inv_w[a] = (-1)*(Pol->wall[i]->word->word[inv_w[0]+1-a]); 
	       if(Pol->wall[t-1]->word->word != NULL)
		 free(Pol->wall[t-1]->word->word);
	       Pol->wall[t-1]->word->word = inv_w;	
	     }
	     if(Pol->wall[i]->word->trans != NULL){
	       A = copy_mat(inv);
	       for(a=0; a < A->rows-1; a++)
		 A->array.SZ[a][A->cols-1] = 0;
	       B = mat_mul(A,Pol->wall[i]->word->trans);
	       for(a=0; a < A->rows-1; a++)
		 B->array.SZ[a][B->cols-1] = (-1)*(B->array.SZ[a][B->cols-1]);
	       free_mat(Pol->wall[t-1]->word->trans); Pol->wall[t-1]->word->trans = NULL;
	       Pol->wall[t-1]->word->trans = B;
	       free_mat(A); A = NULL;
	     }/* if(Pol->wall[t-1]->word->trans != NULL) */
	   }/* if(mat_comp(inv,Pol->wall[t-1]->mat) != 0) */
	 }/* if(i != t-1) */
       }/* if(isok == -1) */
       else{
	 /*   check[t-1] = 0; */
	 old_P_wall_no = Pol->wall_no;
	 incase = 0;
	 unterscheide(i,t-1,Pol,w_v_anz,w_vert_nr,inv,&incase);
	 if(incase == 1){
	   free(w_v_anz);
	   for(k=0; k<Pol->wall_no; k++)
	     free(w_vert_nr[k]);
	   free(w_vert_nr);
	   
	   divide(Pol,t-1,i,inv); /* unterteile w(i) mit Bild von w(t-1) */
	   
	   facts(Pol,&w_v_anz,&w_vert_nr);
	   Bary = bary_list(Pol);
	   
	   check = (int*)realloc(check, Pol->wall_no * sizeof(int));
	   check[Pol->wall_no-1] = 0;
	 }
	 if (incase == 2){
	   free(w_v_anz);
	   for(k=0; k<Pol->wall_no; k++)
	     free(w_vert_nr[k]);
	   free(w_vert_nr);
	   
	   divide(Pol,i,t-1,Pol->wall[i]->mat); /* unterteile w(t-1) mit Bild von w(i) */
	   
	   facts(Pol,&w_v_anz,&w_vert_nr);
	   Bary = bary_list(Pol);
	   
	   check = (int*)realloc(check, Pol->wall_no * sizeof(int));
	   check[Pol->wall_no-1] = 0;
	 }
	 if(old_P_wall_no == Pol->wall_no){
	   printf(" t = %d, i = %d incase = %d \n",t,i,incase);
	   check[i] = -1;
	 }
	 if (incase == 0){
	   check[i] = 1;
	 }
       } /* isok != -1 */
     } /* tmp auf Rand von Pol  und zwar auf Wand t-1*/
     
     if(t == 0){
       old_P_wall_no = Pol->wall_no;
       s = look_for_wall(Pol,i,&zero,&z,w_vert_nr, w_v_anz, &incase); 
       if(incase == 1){
	 unterscheide(i,s,Pol,w_v_anz,w_vert_nr,inv,&incase);
	 if(incase == 1){
	   free(w_v_anz);
	   for(k=0; k<Pol->wall_no; k++)
	     free(w_vert_nr[k]);
	   free(w_vert_nr);
	   
	   divide(Pol,s,i,inv); /* unterteile w(i) mit Bild von w(s) */
	   
	   facts(Pol,&w_v_anz,&w_vert_nr);
	   Bary = bary_list(Pol);
	   
	   check = (int*)realloc(check, Pol->wall_no * sizeof(int));
	   check[Pol->wall_no-1] = 0;
	 }
	 if (incase == 2){
	   free(w_v_anz);
	   for(k=0; k<Pol->wall_no; k++)
	     free(w_vert_nr[k]);
	   free(w_vert_nr);
	   
	   divide(Pol,i,s,Pol->wall[i]->mat); /* unterteile w(s) mit Bild von w(i) */
	   
	   facts(Pol,&w_v_anz,&w_vert_nr);
	   Bary = bary_list(Pol);
	   
	   check = (int*)realloc(check, Pol->wall_no * sizeof(int));
	   check[Pol->wall_no-1] = 0;
	 }
	 if(old_P_wall_no == Pol->wall_no){
	   check[i] = -1;
	 }
	 if (incase == 0){
	   check[i] = 1;
	 }
       }/** if(look_for_wall -> incase == 1) **/
       else {
	 isok = check_vertices(Pol,i,s,w_vert_nr, w_v_anz,inv);
	 if(isok == -1){
	   check[i] = 1; 
	   if(i != t-1){
	     check[t-1] = 1;
	     free_mat(Pol->wall[t-1]->mat);
	     Pol->wall[t-1]->mat = copy_mat(inv);
	   }
	 }
	 else{
	   printf("Fehler in Check_side_trafo \n"); 
	   exit(3);
	 }
       }
     } /* t == 0 */
     
     free_mat(inv); inv = NULL;
     free_mat(tmp); tmp = NULL;
     
     for(i=0; i<Pol->wall_no; i++)
       if(check[i] == 0)
	 break;
     free(zero);
     free(bar);
   }while(i < Pol->wall_no);
   
   
   for(i=0; i<Pol->wall_no; i++){
     check[i] ++; 
     if(check[i] == 10)
       {printf("Fehler \n"); exit(3); }
   }
   
   for(i=0; i<Pol->wall_no; i++)
     if(check[i] == 0)
       break;
   
 }while(i < Pol->wall_no);
 
 free(w_v_anz);
 for(k=0; k<Pol->wall_no; k++)
   free(w_vert_nr[k]);
 free(w_vert_nr);
 
 for(i=0; i<Pol->wall_no; i++)
   free(Bary[i]);
 free(Bary);
 free(check); 
 free_mat(test); test = NULL;
}/** check_side_trafo() **/

static void Trans_anwenden(Pol,length,orb,G,G_inv,Form,P_anz,M, P_inv, P_list)
polyeder_TYP *Pol;
int length;
ele_TYP **orb;
bravais_TYP *G;
matrix_TYP **G_inv;
matrix_TYP *Form;
int P_anz;
matrix_TYP *M;
matrix_TYP **P_inv;
matrix_TYP **P_list;
{
 int		i,j,s,t;
 int 		anz, weiter, neu_anz;
 int 		*rel_wand;
 ele_TYP 	**kand;
 ele_TYP 	**neu_kand;

/*
 P_inv  = (matrix_TYP **)malloc(Pol->wall_no * sizeof(matrix_TYP*));
 P_list  = (matrix_TYP **)malloc(Pol->wall_no * sizeof(matrix_TYP*));
 for(j=0; j<Pol->wall_no; j++){
    P_inv[j] = mat_inv(Pol->wall[j]->mat);	
    P_list[j] = copy_mat(Pol->wall[j]->mat);
 }
 P_anz = Pol->wall_no;
*/
 anz = 1 + (length-1) * P_anz;

 if(anz-1 != 0){
   kand = (ele_TYP **)malloc(anz * sizeof(ele_TYP*));
   for(i=0; i<anz; i++){
     kand[i] = (ele_TYP *)malloc(sizeof(ele_TYP));
   }
   kand[0]->elm = copy_mat(orb[0]->elm);
   kand[0]->trans = NULL;
   kand[0]->schreier_vec = (int*)malloc(sizeof(int));
   kand[0]->schreier_vec[0] = 0;
   s = 1;
   for(i=1; i<length; i++){
     for(j=0; j<P_anz; j++){
       kand[s]->elm = mat_mul(P_list[j],orb[i]->elm);
       Check_mat(kand[s]->elm);
       kand[s]->schreier_vec = 
	 (int*)malloc((orb[i]->schreier_vec[0]+1) * sizeof(int));
       memcpy(kand[s]->schreier_vec, orb[i]->schreier_vec, 
	      (orb[i]->schreier_vec[0]+1)*sizeof(int));
       kand[s]->trans = mat_mul(P_list[j],orb[i]->trans);
       Check_mat(kand[s]->trans);
       s++;
     }
   }/* for(i) */
 }/*if(anz-1!=0)*/
 
/*berechne waende die durch orb und kand geliefert werden und schneide ab
  Polyeder aendert sich !!! */
  
 weiter = 0;
 if(anz-1 != 0){
   for(t=0; t<Pol->wall_no; t++)
     Pol->wall[t]->neu = 1;
   
   weiter = wand(G->gen_no, G->gen, G_inv, kand, anz, Pol, Form, &rel_wand);
   
   while(weiter != 0){
     weiter = 0;
     kand_sort(kand,rel_wand,&anz);
     free(rel_wand);
     neu_kand = create_new_cand(kand,anz,P_list,P_inv,Pol,P_anz); 
     
     for(i=0; i<anz; i++){
       if(kand[i]->elm != NULL){free_mat(kand[i]->elm); kand[i]->elm = NULL;}
       if(kand[i]->trans != NULL){free_mat(kand[i]->trans); kand[i]->trans = NULL;}
       if(kand[i]->schreier_vec != NULL){free(kand[i]->schreier_vec);}
       free(kand[i]); 
     }
     free(kand);
     
     neu_anz = 1 + P_anz * (anz-1) ;
     for(t=0; t<Pol->wall_no; t++)
       Pol->wall[t]->neu = 1;
     
     weiter = 
       wand(G->gen_no, G->gen, G_inv, neu_kand, neu_anz, Pol,Form,&rel_wand);
     
     kand = neu_kand;
     anz = neu_anz;
   }/*while () */
   
   free(rel_wand);
   for(i=0; i<anz; i++){
     if(kand[i]->elm != NULL){free_mat(kand[i]->elm); kand[i]->elm = NULL;}
     if(kand[i]->trans != NULL){free_mat(kand[i]->trans);kand[i]->trans = NULL;}
     if(kand[i]->schreier_vec != NULL){free(kand[i]->schreier_vec);} 
     free(kand[i]);
   }
   free(kand);
   
 }/*if(anz!=0)*/
 /*
   for(i=0; i<P_anz; i++){
   free_mat(P_inv[i]); P_inv[i] = NULL;
   free_mat(P_list[i]); P_list[i] = NULL;
   }
   free(P_inv); free(P_list);
 */
}/** Trans_anwenden () **/

/* Eingabe: G enthaelt zusaetzlich die Erzeuger des Translationsnormalteilers
   als Gruppenerzeuger und zwar sind dies die ersten dim Erzeuger
   in geordneter Reihenfolge (e1, e2, e3 ...,ed).
*/
/**************************************************************************\
@--------------------------------------------------------------------------
@ polyeder_TYP *fub(matrix_TYP *M, bravais_TYP *G, matrix_TYP *F)
@
@ matrix_TYP *M   affine vector, starting point of the algorithm, i.e. (0,0,0,1)
@ bravais_TYP *G  generators of a spacegroup G with translations equal to Z^n,
@                 the first dim generators are the (lexicographically ordered) 
@                 translations.
@ matrix_TYP *F   a positive definite, quadratic form, invariant under the 
@                 pointgroup of G
@
@ calculates a fundamental polyeder Pol for the spacegroup G, the side-
@ transformations for Pol generate G and can be used to determine a 
@ presentation of G
@
@--------------------------------------------------------------------------
\**************************************************************************/
polyeder_TYP *fub(matrix_TYP *M, bravais_TYP *G, matrix_TYP *Form)
{
 int		i,j,t;
 int		kk,k;
 int		P_anz;
 int		dim;
 int		length,leng;
 int		dummy;
 int		weiter;
 int		*rel_wand;
 ele_TYP	**orb;
 ele_TYP	**sorbit;
 matrix_TYP 	*Y,*B;
 matrix_TYP	*I;
 matrix_TYP	**G_inv;
 matrix_TYP	**P_inv, **P_list;
 bravais_TYP	*Stab, *neu_Stab;
 Stab_word_TYP 	*SW,*neu_SW;
 polyeder_TYP 	*Pol;
 polyeder_TYP 	*erg;

extern	ele_TYP **easy_torus_orbit();
extern	ele_TYP **stab_orbit();
extern  void check_side_trafo();

 Stab = init_bravais(G->dim);
 
 SW = NULL;
 
 dim = G->dim;
 I = init_mat(G->dim-1, G->dim-1,"1");
 
 Pol = T_polyeder(M, I, Form);		/* Pol = D(x,Tx) mit Seitentrafo*/
 
 j = check_polyeder(Pol);
 if(j != 0){
   printf("T_poly falsch \n");
   exit(3);}
 
 P_inv  = (matrix_TYP **)malloc(Pol->wall_no * sizeof(matrix_TYP*));
 P_list  = (matrix_TYP **)malloc(Pol->wall_no * sizeof(matrix_TYP*));
 for(j=0; j<Pol->wall_no; j++){
   P_inv[j] = mat_inv(Pol->wall[j]->mat);	
   P_list[j] = copy_mat(Pol->wall[j]->mat);
 }
 P_anz = Pol->wall_no;

 G_inv  = (matrix_TYP **)malloc(G->gen_no * sizeof(matrix_TYP*));
 for(j=0; j<G->gen_no; j++)
   G_inv[j] = mat_inv(G->gen[j]);
 
 /*druecke die Seitentrafo's in den Erzeugern aus (benutzte das es Z^n ist)*/
 for(i=0; i< Pol->wall_no; i++){
   Pol->wall[i]->word = init_word(dim);
   Pol->wall[i]->word->trans = init_mat(dim,1,"");
   for(j=0; j<dim; j++)
     Pol->wall[i]->word->trans->array.SZ[j][0] = 
       Pol->wall[i]->mat->array.SZ[j][dim-1];
 }
 
 orb = easy_torus_orbit( M, G, Stab, &SW, I, &length,Pol,P_inv);


/*transformiere das Bahnelement elm in den Dirichletbereich D(x,Tx) durch 
  sukzessives Anwenden der Seitentransformationen, merke mir dies in trans
 */

 for(i=1; i<length; i++)
   transform(&orb[i]->elm, Pol, &orb[i]->trans, P_inv);

 /*berechne waende die durch orb geliefert werden und schneide ab
   Polyeder aendert sich !!! */

  if(length > 1){
    for(t=0; t<Pol->wall_no; t++)
      Pol->wall[t]->neu = 1;
    
    weiter = wand(G->gen_no, G->gen, G_inv, orb, length, Pol, Form, &rel_wand);
    
    free(rel_wand);
  }

 /* zu jedem orb->elem addiere noch die Translationen der Seiten von DT 
    d.h. wende die entsprechende Seitentrafo an  */

  Trans_anwenden(Pol,length,orb,G,G_inv,Form,P_anz,M,P_inv,P_list);
 
 /* Freigeben der nicht mehr benoetigten Parameter */

  check_side_trafo(Pol,G->gen,G->gen_no,G_inv,M);

  for(i=0; i<P_anz; i++){
    free_mat(P_inv[i]); P_inv[i] = NULL;
    free_mat(P_list[i]); P_list[i] = NULL;
  }
  free(P_inv); free(P_list);
  for(i=0; i<length; i++)
    free_ele(orb[i]);
  free(orb);
  
  B = copy_mat(M);
  while(Stab->gen_no != 0){
    
/* create new element y != x in P . bilde den Mittelpunkt der Verbindungs-
   linie zwischen x und einem Eckpunkt v => y = x+(1/2)*(v-x) = 1/2*(x+v) */

    /* the function random is not really random */
    kk = random_own();

    k = kk % Pol->vert_no;
    Y = mkpt(B,Pol->vert[k]->v);
    
/*bestimme Bahn von y unter Stab, sowie Stab(y) und zugeh. Worte */

    neu_Stab = init_bravais(G->dim);
    neu_SW = NULL;
    leng = 0;

    sorbit = stab_orbit(Y,Stab,SW,neu_Stab,&neu_SW,&leng,G->gen,G_inv,G->gen_no); 

    for(t=0; t<Pol->wall_no; t++)
      Pol->wall[t]->neu = 1;
    
    /*berechne neue Waende und schneide ab*/
    
    dummy = wand(G->gen_no, G->gen, G_inv, sorbit, leng, Pol, Form, &rel_wand);
    
    free(rel_wand);

    check_side_trafo(Pol,G->gen,G->gen_no,G_inv,Y);
    
/* Freigeben der nicht mehr benoetigten Parameter */
    for(i=0; i<Stab->gen_no; i++){
      free(SW[i].sword);
      free_mat(SW[i].trans); SW[i].trans = NULL;
    }
    free(SW); SW = NULL;
    free_bravais(Stab);
    Stab = neu_Stab;
    SW = neu_SW;
    for(i=0; i<leng; i++)
      free_ele(sorbit[i]);
    free(sorbit);
    
    free_mat(B); B = NULL;
    B = Y;
  }/** while(Stab->gen_no != 0) **/
  free_mat(B); B = NULL;
  
  for(i=0; i<Stab->gen_no; i++){
    free(SW[i].sword);
    free_mat(SW[i].trans); SW[i].trans = NULL;
  }
  if(SW != NULL)
    free(SW); SW = NULL;
  free_bravais(Stab); Stab = NULL;
  
  for(i=0; i<G->gen_no; i++){
    free_mat(G_inv[i]); G_inv[i] = NULL;
  }
  free(G_inv);
  
  free_mat(I); I = NULL;
  
  erg = Pol; 
  return(erg);
}/** fub() **/
