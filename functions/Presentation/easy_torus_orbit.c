
/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: easy_torus_orbit.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

#include "typedef.h"
#include "matrix.h"
#include "symm.h"
#include "getput.h"
#include "bravais.h"
#include "sort.h"
#include "presentation.h"

static matrix_TYP *affine_vector_add( matrix_TYP *A, matrix_TYP *B, int	n);
static int isinpoly(matrix_TYP *Vec, polyeder_TYP *P);
static matrix_TYP *transform( matrix_TYP *elm, matrix_TYP **trans, polyeder_TYP *P, matrix_TYP **P_inv);
static int vec_vergleich( matrix_TYP *A, matrix_TYP *B, polyeder_TYP *P, matrix_TYP **P_inv);
static int is_neu( matrix_TYP *A, ele_TYP **orbit, int count, polyeder_TYP *P, matrix_TYP **P_inv);
static void matrix_speichern( matrix_TYP *mat, matrix_TYP ***L, int *listsize,int *anz);
static int mat_vergleich( matrix_TYP *A, matrix_TYP *B);
static struct baum *addbaum( matrix_TYP *mat, matrix_TYP **L, int anz, struct baum *verz, int *schonda);
static matrix_TYP *affine_vector_diff( matrix_TYP *A, matrix_TYP *B);
static int append_new_element(ele_TYP **orbit, matrix_TYP *trans,int orb_length, matrix_TYP *test, int node, int gen_no, matrix_TYP *Erz);


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

static int isinpoly(matrix_TYP *Vec, polyeder_TYP *P)
{
  int 	i, j, s;
  matrix_TYP 	*test,
    *M;
  
  j = Vec->rows;
  s = 0;
  
  M = init_mat(P->wall_no, j, "");
  
  for(i=0; i<M->rows; i++)
    memcpy(M->array.SZ[i], P->wall[i]->gl, j* sizeof(int));
  
  test = mat_mul(M,Vec);
  Check_mat(test);
  
  for (i=0; i<test->rows ; i++){
    if(test->array.SZ[i][0] < 0){
      free_mat(test); test = NULL;
      free_mat(M); M=NULL;
      return(i);	       	/* i = Nr der Wandgl. die verletzt ist*/
    }
  }
  free_mat(test); test = NULL;
  free_mat(M); M=NULL;
  return(-1); 				/*ist im Innern oder auf dem Rand*/ 
}/** isinpoly() **/

/*transformiere das Bahnelement elm in den Dirichletbereich D(x,Tx) durch 
  sukzessives Anwenden der Seitentransformationen, merke mir dies in trans */
static matrix_TYP *transform(elm, trans, P, P_inv)
matrix_TYP *elm;
matrix_TYP **trans;
polyeder_TYP *P;
matrix_TYP **P_inv;
{
  int		test;
  matrix_TYP	*tmp,*tmp1, *tmp2;
  matrix_TYP	*erg;
  
  erg = copy_mat(elm);
  test = isinpoly(erg, P);
  while(test != -1){
    tmp = mat_mul(P_inv[test], erg);
    Check_mat(tmp);
    free_mat(erg); erg = NULL;
    erg = tmp;
    test = isinpoly(erg, P);
  }
  tmp1 = affine_vector_add(erg,elm,-1);  /* gx+t = g'x = erg */
  
  if(trans[0] != NULL){
    tmp2 = affine_vector_add(trans[0],tmp1,1); /*addiere tmp1 zu trans hinzu*/
    free_mat(trans[0]); trans[0] = NULL;
    trans[0] = tmp2; 
    free_mat(tmp1), tmp1 = NULL;
  }
  else
    trans[0] = tmp1; 
  return(erg);
}/** transform() **/


static int vec_vergleich(A, B,P,P_inv)
matrix_TYP *A, *B ;
polyeder_TYP *P;
matrix_TYP **P_inv;
{
  int s;

  s = mat_comp(A,B);  
  if(s == 0){
    return(0);
  }
  else{
    return(1);
  }
  
}/**vec_vergleich()**/

static int is_neu(A, orbit, count,P,P_inv)
matrix_TYP *A;
ele_TYP **orbit;
int count;
polyeder_TYP *P;
matrix_TYP **P_inv;
{
  int i;
  int schonda;
  schonda = 1;
  i = 0;
  
  while(schonda == 1 && i < count) {
    schonda = vec_vergleich(A, orbit[i]->elm,P,P_inv);
    i++;
  }
  if(schonda == 1)
    schonda = -1;
  if(schonda == 0)	/*war schonmal da*/
    schonda = i-1;
  return(schonda);
}/**is_neu()**/

static void matrix_speichern(mat, L, listsize, anz)
matrix_TYP *mat;
matrix_TYP ***L;
int *listsize, *anz;
{
  if((*anz) == 0)
    if ((*L=(matrix_TYP **)malloc((*listsize)*sizeof(matrix_TYP *)))==NULL)
      {
	fprintf (stderr, "malloc failed\n");
	exit (2);
      }
  if ((*anz) >= (*listsize))
    {
      *listsize += EXT_SIZE;
      if ((*L=(matrix_TYP **)realloc(*L,(*listsize)*sizeof(matrix_TYP *)))==NULL)
	{
	  fprintf (stderr, "realloc failed\n");
	  exit (2);
	}
    }
  (*L)[(*anz)] = mat; 
  (*anz) ++;
}/**matrix_speichern()**/

static int mat_vergleich(A, B)
matrix_TYP *A, *B;
{
  int i,j;
  
  for(i=0; i<A->rows; i++){
    for(j=0; j<A->cols; j++)
      if(B->kgv * A->array.SZ[i][j] != A->kgv * B->array.SZ[i][j])
	{
          if(B->kgv * A->array.SZ[i][j] < A->kgv * B->array.SZ[i][j])
	    return(-1);
          return(1);
	}
  } 
  return(0);
}/**mat_vergleich()**/

static struct baum *addbaum(mat, L, anz, verz, schonda)
matrix_TYP *mat;
matrix_TYP **L;
int anz;
struct baum *verz;
int *schonda;
{
  int vergleich;
  *schonda = -1;
  if(verz == NULL)
    {
      verz = (struct baum *) malloc(sizeof(struct baum));
      verz->no = anz;
      verz->left = verz->right = NULL;
    }
  else
    {
      vergleich = mat_vergleich(mat, L[verz->no]);
      if(vergleich<0)
	verz->left = addbaum(mat, L, anz, verz->left, schonda);
      if(vergleich > 0)
	verz->right = addbaum(mat, L, anz, verz->right, schonda);
      if(vergleich == 0)
	*schonda = verz->no;
    }
  return(verz);
}/**addbaum()**/

static matrix_TYP *affine_vector_diff(A,B) /* SPALTEN-KONVENTION*/
matrix_TYP *A;
matrix_TYP *B;
{  
  matrix_TYP *C;
  int i,j;

  C = init_mat(A->rows, A->cols,"");
  for(i=0; i<C->rows-1; i++){
    for(j=0; j<C->cols; j++){
      C->array.SZ[i][j] = 
	B->kgv*A->array.SZ[i][j]-A->kgv*B->array.SZ[i][j];
    }
  }
  C->kgv = A->kgv * B->kgv;
  i = C->rows; 
  j = C->cols;
  C->array.SZ[i-1][j-1] = C->kgv;
  Check_mat(C);
  return(C);
}/***affine_vector_diff()***/

static int append_new_element(ele_TYP **orbit, matrix_TYP *trans,int orb_length,
			      matrix_TYP *test, int node, int gen_no, matrix_TYP *Erz)
{
  int        	i;
  int   	l,k,SIZE;
  matrix_TYP 	*A, *B;
  
 SIZE = 256;
 	
 if(orb_length == SIZE - 1){
   SIZE = SIZE+256;
   orbit=(ele_TYP  **)realloc(orbit,SIZE *sizeof(ele_TYP *));
   for(k=0; k < 256; k++)
     orbit[SIZE-256+k] = (ele_TYP *)malloc(1 *sizeof(ele_TYP));
 }
 
 if(orbit[node]->trans != NULL)
   {
     A = copy_mat(Erz);
     
     for(i=0; i<A->rows-1; i++)
       A->array.SZ[i][A->cols-1] = 0;
     B = mat_mul(A,orbit[node]->trans);
     Check_mat(B);
     free_mat(A); A = NULL;
     orbit[orb_length]->trans = affine_vector_add(trans,B,1);
     free_mat(B); B = NULL; 
   }
 else
   {
     orbit[orb_length]->trans = copy_mat(trans);
   }
 l = orbit[node]->schreier_vec[0]+1;
 orbit[orb_length]->elm = test;
 orbit[orb_length]->schreier_vec = (int*)calloc((l+1),sizeof(int));
 
 orbit[orb_length]->schreier_vec[0]= orbit[node]->schreier_vec[0]+1;
 
 if(l-1 != 0){
   for(i=1;i<l;i++)
     orbit[orb_length]->schreier_vec[i+1] = orbit[node]->schreier_vec[i];
 }
 orbit[orb_length]->schreier_vec[1] = gen_no+1; 
 /* Index +1, spaeter beachten*/
 return(SIZE);
}/* function append_new_element */

/*erzeuge Wort, welches das Stabilisator-Element in den Erz der Gr. ausdr.*/
static void Word(erz_no,orb_no,SW, Stab_gen_no, orbit, schonda, B2)
int 		erz_no,
    		orb_no;
Stab_word_TYP   **SW;
int 		Stab_gen_no;
ele_TYP 	**orbit;
int 		schonda;
matrix_TYP 	*B2;
{
  int 	a,b,c;
  int 	i,j;
  int   s;
  
  j = erz_no;
  i = orb_no;

  if(*SW == NULL){
    *SW = (Stab_word_TYP *)malloc(Stab_gen_no * sizeof(Stab_word_TYP));
    for(s=0; s < Stab_gen_no; s++)
      (*SW)[s].trans = (matrix_TYP *)malloc(sizeof(matrix_TYP));
  }
  else{
    *SW = 
      (Stab_word_TYP *)realloc((*SW),Stab_gen_no * sizeof(Stab_word_TYP));
    (*SW)[Stab_gen_no-1].trans = (matrix_TYP *)malloc(sizeof(matrix_TYP));
  }
  
  a = orbit[schonda]->schreier_vec[0]+orbit[i]->schreier_vec[0]+1;
  (*SW)[Stab_gen_no-1].sword = (int*)malloc((a+1)*sizeof(int));
  
  a = orbit[schonda]->schreier_vec[0]; 
  for(b=0,c=1; b<a; b++){
    (*SW)[Stab_gen_no-1].sword[c] = 
      (-1)*(orbit[schonda]->schreier_vec[a-b]); 
    c++;
  }
  (*SW)[Stab_gen_no-1].sword[c++] = j+1;
  
  a = orbit[i]->schreier_vec[0];
  for(b=0; b<a; b++){
    (*SW)[Stab_gen_no-1].sword[c++] =
      orbit[i]->schreier_vec[b+1]; 
  }
  (*SW)[Stab_gen_no-1].sword[0] = c-1;
  
  (*SW)[Stab_gen_no-1].trans = B2;
}/**Word**/

/**************************************************************************\
@--------------------------------------------------------------------------
@ ele_TYP **easy_torus_orbit( M, G, Stab,SW,Tmat, l, P, P_inv)
@
@ matrix_TYP *M;
@ bravais_TYP *G;	 given spacegroup G
@ bravais_TYP *Stab;     will be the stabilizer of M under G 
@ Stab_word_TYP   **SW;  SW = NULL, will contain the words of Stab expressed in @                        the generators of G
@ matrix_TYP *Tmat;      grammatrix of the translationlattice
@ int *l;                l will denote the length of the orbit
@ polyeder_TYP *P;       P denotes a polyeder
@ matrix_TYP **P_inv;    a list of the inverse sidetransformations of P
@
@ calculates the orbit of M under G in the polyeder P and the stabilizer of M
@ as well as the words for expressing Stab in the generators for G 
@--------------------------------------------------------------------------
\**************************************************************************/
ele_TYP **easy_torus_orbit( M, G, Stab,SW,Tmat, l, P, P_inv)
matrix_TYP *M;
bravais_TYP *G;
bravais_TYP *Stab;
Stab_word_TYP   **SW;
matrix_TYP *Tmat;
int *l;
polyeder_TYP *P;
matrix_TYP **P_inv;
{
  ele_TYP 		**orbit;
  int			i,j;
  int 			SIZE;
  int 			count,old;
  int 			neu_anz;
  int 			schonda, sch;
  int 			histanz, histinvanz, histsize, histinvsize;
  int 			erz_anz, Ssize;
  int 			size;
  rational		re,li;
  matrix_TYP 		*stabtest;
  matrix_TYP 		*test,*test1;
  matrix_TYP 		*I;
  matrix_TYP 		*I1, *I2;
  matrix_TYP   		*K1, *K2, *Ih, *K2h;
  matrix_TYP 		**hist, **histinv;
  matrix_TYP 		**Erz, **Erz_inv;
  matrix_TYP		*trans;
  struct baum 		*Sverz;
  matrix_TYP 	*A, *B1, *sst, *STS; 
  int 		u;

  SIZE = 256;
  orbit = (ele_TYP **)malloc(SIZE *sizeof(ele_TYP*));
  
  orbit[0] = (ele_TYP *)malloc(1 *sizeof(ele_TYP));
  orbit[0]->schreier_vec = (int*)malloc(1*sizeof(int));
  orbit[0]->schreier_vec[0] = 0;
  orbit[0]->trans = NULL;
  
  for(i=1; i<SIZE; i++){
    orbit[i] = (ele_TYP *)malloc(1 *sizeof(ele_TYP));
    orbit[i]->schreier_vec = NULL;
    orbit[i]->trans = NULL;
  }
  /*schreier_vec[0] = Laenge des Wortes*/
  erz_anz = G->gen_no;
  Ssize = 0;
  size = 0;
  histsize = 0;
  histinvsize = 0;
  histanz = 0;
  histinvanz = 0;
  *l = 1;
  old = 0;
  count = *l;
  Sverz = NULL;
  Erz = (matrix_TYP **) malloc(erz_anz *sizeof(matrix_TYP *));
  for(i=0;i< G->gen_no;i++)
    Erz[i] = G->gen[i];
  
  Erz_inv = (matrix_TYP **) malloc(erz_anz *sizeof(matrix_TYP *));
  for(i=0;i<erz_anz;i++)
    Erz_inv[i] = mat_inv(Erz[i]);
  
  I = einheitsmatrix(G->gen[0]->rows);
  I1 = einheitsmatrix(G->gen[0]->rows);
  I2 = einheitsmatrix(G->gen[0]->rows);
  matrix_speichern(I1,&hist,&histsize,&histanz);
  matrix_speichern(I2,&histinv,&histinvsize,&histinvanz);
  trans = NULL;
  orbit[0]->elm = copy_mat(M);

  re.z = -1; re.n = 1; li.z = 1; li.n = 1;
  neu_anz = 1;
  while(neu_anz !=0){
    neu_anz = 0;
    
    for(i=old; i<*l; i++){
      for(j=0; j<erz_anz; j++){ 
	test1= mat_mul( Erz[j], orbit[i]->elm);
	Check_mat(test1);
	test = transform( test1, &trans, P, P_inv);
	Check_mat(test);
	free_mat(test1); test1 = NULL;
	schonda = is_neu( test, orbit, count, P, P_inv); 
	if(schonda == -1){
	  size = append_new_element( orbit, trans, count, test, i, j, Erz[j]);
          
	  free_mat(trans); trans = NULL;
	  K1 = mat_mul(Erz[j], hist[i]); 	/*anne, 15.8. (SPALTEN)*/
	  matrix_speichern(K1, &hist, &histsize, &histanz);
	  K1 = NULL;
	  K1 = mat_mul(histinv[i],Erz_inv[j]);
	  matrix_speichern(K1, &histinv, &histinvsize, &histinvanz);
	  K1 = NULL;
	  
	  count++;
	  neu_anz++;
	}/*schonda*/
	else{ 
	  K1 = mat_mul(Erz[j], hist[i]);
	  K2 = mat_mul(histinv[schonda], K1);
	  K2h = copy_mat(K2); 
	  real_mat(K2h,K2h->rows-1,K2h->cols-1);
	  Ih = einheitsmatrix(K2h->rows); 
	  
	  if(mat_vergleich(K2h, Ih) != 0)
	    {
	      /* berechne tatsaechlichen Stabilisator vom Punkt M */
	      
	      A = mat_mul(K2,M);
	      if(mat_vergleich(A,M) != 0){
		B1 = affine_vector_diff(M,A); 
		free_mat(A); A = NULL;
		sst = init_mat(K2->rows,K2->cols,"1");
		sscal_mul(sst,B1->kgv);
		sst->kgv = B1->kgv;
		for(u=0; u<sst->rows-1; u++)
		  sst->array.SZ[u][sst->cols-1] = B1->array.SZ[u][0];
		Check_mat(sst);
		STS = mat_mul(sst,K2); 
		Check_mat(STS);
		free_mat(sst); sst = NULL;
		
		/* merke mir B1 beim Stabword 16.6.97*/

		Sverz = addbaum(STS, Stab->gen, Stab->gen_no, Sverz, &sch);
		if(sch == -1){
		  free_mat(K2); K2 = NULL;
		  matrix_speichern(STS, &Stab->gen, &Ssize, &Stab->gen_no);
		  STS = NULL; 

		  Word(j,i,SW, Stab->gen_no, orbit, schonda, B1);
		  stabtest = Stab_word_to_mat((*SW)[Stab->gen_no-1], Erz, Erz_inv, erz_anz,Tmat);
		  free_mat(stabtest); stabtest = NULL;
		  
		}/***if(sch == -1)***/
	      }/*if(mat_vergleich(A,M) != 0)*/
	      else{
                Sverz = addbaum(K2, Stab->gen, Stab->gen_no, Sverz, &sch);
              	if(sch == -1){
		  matrix_speichern(K2, &Stab->gen, &Ssize, &Stab->gen_no);
		  
		  B1 = init_mat(K2->rows, 1, "");
		  B1->array.SZ[B1->rows-1][0] = 1; 
		  Word(j,i,SW, Stab->gen_no, orbit, schonda, B1);
		  stabtest = Stab_word_to_mat((*SW)[Stab->gen_no-1], Erz, Erz_inv, erz_anz,Tmat);
		  free_mat(stabtest); stabtest = NULL;
		}/**if(sch == -1) **/
	       }/*else*/
   	      
	    }/***if(mat_vergleich(K2h, Ih) != 0)***/
	  
	  if(schonda == -1 || (schonda != -1 && sch != -1)){
	    free_mat(K2);  K2 = NULL;
	  }  
	  free_mat(K1); K1 = NULL; K2 = NULL;
	  free_mat(Ih); Ih = NULL;
	  free_mat(K2h); K2h = NULL;
	  
	  if(trans != NULL){
	    free_mat(trans); trans = NULL;
	  }
	  free_mat(test); test = NULL;
	}/**else **/
      }/*for(j)*/
    }/*for(i)*/
    old = *l;
    *l = count;
  }/*while*/
  
  for(i=0; i<histanz; i++)
    {
      free_mat(hist[i]); hist[i] = NULL;
      free_mat(histinv[i]); histinv[i] = NULL;
    }
  free(hist); free(histinv);
  free_mat(I); I = NULL;
  
  for(i=0;i<erz_anz;i++){
    free_mat(Erz_inv[i]); Erz_inv[i] = NULL;
  }
  free(Erz); free(Erz_inv);
  
  for(i=(*l); i<size; i++){
    free(orbit[i]);
  }
  
  orbit=(ele_TYP  **)realloc(orbit,(*l) *sizeof(ele_TYP *));
  
  return(orbit);
}/* easy_torus_orbit() */
