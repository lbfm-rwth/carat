
/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: Back.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

#include"typedef.h"
#include"matrix.h"
#include"polyeder.h"
#include"sort.h"
#include"presentation.h"
#include"getput.h"
#include"symm.h"
#include"bravais.h"

static matrix_TYP *affine_vector_add( matrix_TYP *A, matrix_TYP *B, int n);
static int is_in_Poly(matrix_TYP *test, polyeder_TYP *P);
static int distance(matrix_TYP *A, matrix_TYP *B, matrix_TYP *Form);
static int  search_in_erglist(int a,int *erglist);

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

static int is_in_Poly(matrix_TYP *test, polyeder_TYP *P)
{
  int i,j,s;
  matrix_TYP *M, *tmp;
  
  j = test->rows;
  M = init_mat(P->wall_no,j,""); 
  
  for(i=0; i<M->rows; i++)
    memcpy(M->array.SZ[i], P->wall[i]->gl, j*sizeof(int*));
  
  tmp = mat_mul(M,test);
  Check_mat(tmp);
  
  for(i=0; i<P->wall_no; i++){
    if(tmp->array.SZ[i][0] < 0)
      return(i);
  } 
  s = -1;
  return(s);
} /** is_in_Poly() **/

static int distance(matrix_TYP *A, matrix_TYP *B, matrix_TYP *Form)
{  
  int	d;
  matrix_TYP *tmp,*tmp1, *tmp2,*tmp3;
  
  tmp = affine_vector_add(A,B,-1); Check_mat(tmp);
  tmp->rows --;
  tmp1 = tr_pose(tmp); Check_mat(tmp1);
  tmp2 =  mat_mul(tmp1, Form); Check_mat(tmp2);
  tmp3 = mat_mul(tmp2,tmp); Check_mat(tmp3);
  
  d = tmp3->array.SZ[0][0] * tmp->kgv;
  
  return(d);
}/*distance()*/

static int  search_in_erglist(int a,int *erglist)
{ 
  int b;  
  
  for(b=0; b<erglist[0]; b++){
    if(erglist[b+1] == abs(a) )
      return(1);
  }
  return(0);
}/*search_in_erglist */

/* erg[i] = wort wie sich der i-te Erzeuger in den Seitentrafo's ausdr.
       ACHTUNG hier zaehle ich von 0 ab.(nicht ab 1 wie bei Stab / Pres) 
       vec ist ein generischer Vektor im Innern des Fundamentalbereiches   
   */

/**************************************************************************\
@
@int **back(matrix_TYP *vec, matrix_TYP **Erz, int erz_anz, polyeder_TYP *Pol,
@	   matrix_TYP *Form)
@
@vec is a vector in the inner part of the polyeder Pol  
@Erz is a list of generators of a spacegroup R and erz_anz the number of 
@    generators
@Pol   is a fundamental polyeder for R
@Form is a positive definite form according to which Pol was calculated
@
@returns a list of words how the generators of the group are expressed in 
@terms of the sidetransformations of the fundamental polyeder Pol.
@
\**************************************************************************/
int **back(matrix_TYP *vec, matrix_TYP **Erz, int erz_anz, polyeder_TYP *Pol,
	   matrix_TYP *Form)
{
  int  		i,j,s;
  int		a;
  int		SIZE;
  int		diff;
  int		count;
  int		**erg;
  int		*erglist;
  int		*Liste, *Liste_inv;
  matrix_TYP 	*test;
  matrix_TYP 	*tmp;
  matrix_TYP 	**next;
  matrix_TYP 	**P_inv;
  
  SIZE = 32;
  
  /*next = Liste der Nachbarn von Vec */
  next = (matrix_TYP **)malloc(Pol->wall_no * sizeof(matrix_TYP*));
  for(i=0; i<Pol->wall_no; i++){
    Check_mat(Pol->wall[i]->mat);
    next[i] = mat_mul(Pol->wall[i]->mat,vec);
    Check_mat(next[i]);
  }
  
  P_inv = (matrix_TYP **)malloc(Pol->wall_no * sizeof(matrix_TYP*));
  for(i=0; i<Pol->wall_no; i++){
    P_inv[i] = mat_inv(Pol->wall[i]->mat);
  }
  
  Liste = generate_Liste(Pol,&Liste_inv);
 /* Liste der Erzeuger (besteht gerade aus den P->wall[i]->mat mit Liste[i]>0)*/
  erglist = (int*)calloc(Pol->wall_no+1, sizeof(int));
  s = 1;
  for(i=0; i<Pol->wall_no; i++){
    erglist[s] = i;
    s++;
  }
  erglist[0] = s-1;
  erglist = (int*)realloc(erglist, s*sizeof(int));
  
  /* erg[i] = wort wie sich der i-te Erzeuger in den Seitentrafo's ausdr.*/
  erg = (int**)malloc(erz_anz * sizeof(int*));
  
  for(i=0; i<erz_anz; i++){
    erg[i] = (int*)calloc(SIZE, sizeof(int));   
    
    test = mat_mul(Erz[i],vec);
    diff = distance(test,vec,Form);
    count = 1;
    while(diff != 0){
      s = is_in_Poly(test,Pol);
      if(s != -1){
	if(count == SIZE)
	  { SIZE = count + SIZE;
	  erg[i] = (int*)realloc(erg[i], SIZE *sizeof(int));
	  }
	erg[i][count++] = s;
	tmp = mat_mul(P_inv[s],test);
	free_mat(test); test = NULL;
	test = tmp;
	diff = distance(test,vec,Form);
      }
      else
        {printf("Fehler in Back bei Erzeuger %d \n", i); exit(3);}
    }      
    erg[i][0] = count-1;   /*anne 18.8.*/
    erg[i] = (int*)realloc(erg[i], count *sizeof(int)); /*anne 18.8.*/
    for(j=0; j<erg[i][0]; j++){
      a = erg[i][j+1];
      s = search_in_erglist(abs(a),erglist);
      
      if(s == 0){ /* nicht in erglist */
	/*ersetze a=erg[i][j+1] durch Erzeuger =(Pol->wall[a]->mat)^(-1) = Liste_int[a] */
	if(a > 0) 
	  erg[i][j+1] = Liste_inv[a];
	else
	  erg[i][j+1] = (-1)*Liste_inv[a];
      }/*if(s == 0)*/
    }/*for(j)*/
    
    /* next 7 lines inserted 1.9.97 to count from 1 */
    for(j=0; j<erg[i][0]; j++){
      a = erg[i][j+1];
      if(a < 0)
	erg[i][j+1] = a-1;
      else
	erg[i][j+1] = a+1;
    }/*for(j)*/
    
  }/** for(i=0; i<erz_anz; i++) **/
  
  return(erg);
}/** back() **/
