
/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: stab_orbit.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

#include <stdlib.h>
#include "typedef.h"
#include "matrix.h"
#include "symm.h"
#include "getput.h"
#include "bravais.h"
#include "sort.h"
#include "base.h"
#include "presentation.h"

typedef struct{
matrix_TYP *elm;
int *schreier_vec;
int *neu_schreier_vec;
matrix_TYP *trans;
} neu_ele_TYP;

static int vec_vergleich( matrix_TYP *A, matrix_TYP *B );
static int is_neu( matrix_TYP *A, neu_ele_TYP **orbit, int count); 
static void A_matrix_speichern( matrix_TYP *mat, matrix_TYP ***L, int *listsize, int *anz);
static void matrix_speichern( matrix_TYP *mat, matrix_TYP ***L, int *listsize, int *anz);
static int mat_vergleich( matrix_TYP *A, matrix_TYP *B);
static void Liste( matrix_TYP *mat, matrix_TYP ***List, int *anz, int *schonda, int *Ssize);
static int *trafo( int *schreier_vec, matrix_TYP **trans, bravais_TYP *Stab, Stab_word_TYP *SW, matrix_TYP **Mat, matrix_TYP **Mat_inv, int matanz);
static int append_new_element( neu_ele_TYP ***orbit, int orb_length, matrix_TYP *test, int node, int gen_no, bravais_TYP *Stab, Stab_word_TYP *SW, matrix_TYP 	**Mat, matrix_TYP **Mat_inv, int matanz,int Size);
static void Word( int erz_no,int orb_no, Stab_word_TYP  **SW, int Stab_gen_no, neu_ele_TYP **orbit, int schonda);
static void change( Stab_word_TYP **neu_SW, Stab_word_TYP *SW, int neu_anz, bravais_TYP *Stab, matrix_TYP **Mat,matrix_TYP **Mat_inv, int matanz);
static int Finde( matrix_TYP *mat, bravais_TYP *brav, int l); 

static int vec_vergleich(A, B)
matrix_TYP *A, *B ;
{
  int s,k;
  s = mat_comp(A,B);
  k=abs(s);
  return(k);
}/**vec_vergleich()**/

static int is_neu(A, orbit, count)
     matrix_TYP *A;
     neu_ele_TYP **orbit;
     int count;
{
  int i;
  int schonda;
  schonda = 1;
  i = 0;
  
  while(schonda == 1 && i < count) {
    schonda = vec_vergleich(A, orbit[i]->elm);
    i++;
  }
  if(schonda == 1)
    schonda = -1;
  if(schonda == 0)	/*war schonmal da*/
    schonda = i-1;
  return(schonda);
}/**is_neu()**/

static void A_matrix_speichern(mat, L, listsize, anz)
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
  (*anz) ++;
  *L=(matrix_TYP **)realloc(*L,(*anz)*sizeof(matrix_TYP *));
  (*L)[(*anz-1)] = mat; 
}/**matrix_speichern()**/

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

static void Liste(mat,List,anz,schonda,Ssize)
matrix_TYP *mat;
matrix_TYP ***List;
int *anz;
int *schonda;
int	*Ssize;
{
  int 	i;
  int 	vergleich;
  
  *schonda = -1;
  
  if(*List == NULL){
    *List = (matrix_TYP **)malloc(sizeof(matrix_TYP *));
    List[0][0] = copy_mat(mat);
    *anz = 1;
  }
  else{
    for(i=0; i < *anz; i++){
      vergleich = mat_vergleich(mat, List[0][i]);
      if(vergleich == 0){
	*schonda = i; 
	return;}
    }
    A_matrix_speichern(mat, List, Ssize, anz);
  }/* else */
}/* Liste() */

/* die Worte in schreier_vec noch umschreiben und bezueglich der Erzeuger der
   Raumgruppe ausdruecken */
static int *trafo(schreier_vec, trans, Stab, SW, Mat, Mat_inv, matanz)
int *schreier_vec;
matrix_TYP **trans;
bravais_TYP *Stab;
Stab_word_TYP *SW; 
matrix_TYP **Mat, **Mat_inv;
int matanz;
{
  int		i,j;
  int		a,b,c;
  int  		*erg;
  Stab_word_TYP *tmp, *tmp1;
  
  b = 0;
  for(i=1; i <= schreier_vec[0]; i++){
    a = abs(schreier_vec[i]) - 1;
    b += SW[a].sword[0];
  }
  
  erg = (int *)malloc((b+1) * sizeof(int));
  
  c = 1;
  for(i=1; i <= schreier_vec[0]; i++){
    a = abs(schreier_vec[i]) - 1;
    for(j=1; j <= SW[a].sword[0]; j++){
      if(schreier_vec[i] > 0)
	erg[c++] = SW[a].sword[j];	       	/*Invertiere bei a<0...*/
      else
	erg[c++] = (-1)*(SW[a].sword[(SW[a].sword[0]+1-j)]);
    }
  } 
  c--;
  if(c != b){
    printf("Fehler in trafo \n");
    exit(3);
  }
  erg[0] = c;
  
  /*einsetzen und "trans" vorbeiziehen */ 
  
  if(schreier_vec[0] != 0){ 
    a = abs(schreier_vec[schreier_vec[0]])-1;
    tmp = (Stab_word_TYP*)malloc(sizeof(Stab_word_TYP));
    tmp1 = (Stab_word_TYP*)malloc(sizeof(Stab_word_TYP));
    tmp1[0].sword = (int*)malloc((SW[a].sword[0]+1) * sizeof(int));
    memcpy(tmp1[0].sword, SW[a].sword, (SW[a].sword[0]+1)*sizeof(int));
    tmp1[0].trans = copy_mat(SW[a].trans);
    
    for(i=schreier_vec[0]-1; i>0; i--){
      b = abs(schreier_vec[i])-1;
      
      Stab_word_mul(&SW[b], tmp1, Mat, Mat_inv, matanz, tmp);
      
      free(tmp1[0].sword);
      if(tmp1[0].trans != NULL){
	free_mat(tmp1[0].trans); tmp1[0].trans = NULL;
      }
      
      tmp1[0].sword = (int*)malloc((tmp[0].sword[0]+1) * sizeof(int));
      memcpy(tmp1[0].sword, tmp[0].sword, (tmp[0].sword[0]+1)*sizeof(int));
      tmp1[0].trans = copy_mat(tmp[0].trans);
    }
    trans[0] = tmp1[0].trans;
    tmp1[0].trans = NULL;
    free_Stab_word(tmp1);
  }
  return(erg);
}/** trafo() **/

static int append_new_element( orbit, orb_length, test, node, gen_no, Stab, SW, Mat, Mat_inv, matanz,SIZE)
     neu_ele_TYP 	***orbit;
     int 		orb_length;
     matrix_TYP 	*test;
     int 		node;
     int 		gen_no;
     bravais_TYP 	*Stab; 
     Stab_word_TYP 	*SW;
     matrix_TYP 	**Mat, **Mat_inv;
     int 		matanz;
     int                SIZE;
{
  int        	i;
  int   	l,k;
  int		*new_schreier_vec;

  if(orb_length == SIZE - 1){
    SIZE = SIZE+MIN_SPEICHER;
    orbit[0]=(neu_ele_TYP  **)realloc(orbit[0],SIZE *sizeof(neu_ele_TYP *));
    for(k=0; k < MIN_SPEICHER; k++)
      orbit[0][SIZE-MIN_SPEICHER+k] = (neu_ele_TYP *)calloc(1 ,sizeof(neu_ele_TYP));
  }
  
  l = orbit[0][node]->schreier_vec[0]+1;
  orbit[0][orb_length]->elm = test;
  orbit[0][orb_length]->trans = NULL;
  orbit[0][orb_length]->schreier_vec = (int*)malloc((l+1)*sizeof(int));
  
  orbit[0][orb_length]->schreier_vec[0]= orbit[0][node]->schreier_vec[0]+1;
  
  if(l-1 != 0){
    for(i=1;i<l;i++)
     orbit[0][orb_length]->schreier_vec[i+1] = orbit[0][node]->schreier_vec[i];
  }
  orbit[0][orb_length]->schreier_vec[1] = gen_no+1; 
  /* Index +1, spaeter beachten*/
  
  /* die Worte in schreier_vec noch umschreiben und bezueglich der Erzeuger der
     Raumgruppe ausdruecken */
  
  new_schreier_vec = 
    trafo(orbit[0][orb_length]->schreier_vec, &(orbit[0][orb_length]->trans), Stab, SW,Mat, Mat_inv, matanz);
  orbit[0][orb_length]->neu_schreier_vec = new_schreier_vec;
  
  return(SIZE);
}/* function append_new_element */

/*erzeuge Wort, welches das Stabilisator-Element in den Erz der Gr. (also des 
  vorhergehenden Stabilisators) ausdrueckt*/
static void Word(erz_no,orb_no,SW, Stab_gen_no, orbit, schonda)
     int 		erz_no, orb_no;
     Stab_word_TYP      **SW;
     int 		Stab_gen_no;
     neu_ele_TYP 	**orbit;
     int 		schonda;
{
  int 	a,b,c;
  int 	i,j,s;
  
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
  
  (*SW)[Stab_gen_no-1].trans = NULL;
}/**Word**/

static void change(neu_SW, SW, neu_anz,Stab, Mat, Mat_inv,matanz)
     Stab_word_TYP **neu_SW, *SW;
     int	      neu_anz;
     bravais_TYP *Stab;
     matrix_TYP **Mat,**Mat_inv;
     int matanz;
{
  int j;
  int *erg;
  
  for(j = 0; j< neu_anz; j++){
    erg = trafo(neu_SW[0][j].sword, &(neu_SW[0][j].trans),Stab,SW,Mat,Mat_inv,matanz);
    if(neu_SW[0][j].sword != NULL)
      free(neu_SW[0][j].sword);
    neu_SW[0][j].sword = (int*)malloc((erg[0]+1)*sizeof(int));
    memcpy(neu_SW[0][j].sword, erg,(erg[0]+1)*sizeof(int));
    free(erg);
  }/* for(j) */
  
}/** change() **/

static int Finde(mat, brav ,l)
     matrix_TYP *mat;
     bravais_TYP *brav;
     int l;
{
  int i;
  
  for(i=l; i<brav->gen_no; i++){
    if(mat_comp(mat,brav->gen[i]) == 0)
      return(i);
  }
  printf("Fehler beim Rueckgaengig machen der Permutation\n");
  return(-1);
}/* Finde() */

/**************************************************************************\
@--------------------------------------------------------------------------
@ ele_TYP **stab_orbit()
@ matrix_TYP 	*Y;
@ bravais_TYP 	*Stab, *neu_Stab;
@ Stab_word_TYP   *SW, **neu_SW;
@ int 		*leng;
@ matrix_TYP 	**Mat, **Mat_inv;
@ int 		matanz;
@ 
@ calculates the orbit of a point Y under a group Stab and its stabilizer 
@ neu_Stab. Stab is in general a finite subgroup of a spacegroup R, SW gives the
@ words how the generating elements of Stab are expressed in terms of a 
@ generating  set of R (given in Mat, the inverses in Mat_inv, their number in 
@ matanz).
@ In leng the length of the orbit is stored.
@ One gets for each element of the stabilizer a word expressing it in the given 
@ generators of the group R and stored in neu_Stab. 
@ Also the schreier_vectors of the orbit elements are calculated.
@ 
@--------------------------------------------------------------------------
\**************************************************************************/
ele_TYP **stab_orbit(Y,Stab, SW, neu_Stab, neu_SW, leng, Mat, Mat_inv, matanz)  
matrix_TYP 	*Y;
bravais_TYP 	*Stab, *neu_Stab;
Stab_word_TYP   *SW, **neu_SW;
int 		*leng;
matrix_TYP 	**Mat, **Mat_inv;
int 		matanz;
{
  neu_ele_TYP 		**orbit;
  ele_TYP 		**erg;
  int			i,j,l,ll;
  int 			SIZE;
  int 			count,old;
  int 			neu_anz;
  int 			schonda, sch;
  int 			histanz, histinvanz, histsize, histinvsize;
  int 			erz_anz, Ssize;
  int 			Size;
  int			red, old_gen_no;
  rational		re,li;
  matrix_TYP 		*test;
  matrix_TYP 		*tmp;
  matrix_TYP 		*I;
  matrix_TYP 		*I1, *I2;
  matrix_TYP   		*K1, *K2;
  matrix_TYP 		**base;
  matrix_TYP 		**hist, **histinv;
  matrix_TYP 		**Erz, **Erz_inv;
  matrix_TYP		**merk_stab;
  struct baum 		*Sverz;
  bahn			**strong;

SIZE = MIN_SPEICHER;
	orbit = (neu_ele_TYP **)malloc(SIZE *sizeof(neu_ele_TYP*));
	   for(i=0; i<SIZE; i++){
	      orbit[i] = (neu_ele_TYP *)malloc(1 *sizeof(neu_ele_TYP));
	      orbit[i]->schreier_vec = (int*)malloc(1*sizeof(int));
	      orbit[i]->schreier_vec[0] = 0;
	      orbit[i]->neu_schreier_vec = (int*)malloc(1*sizeof(int));
	      orbit[i]->neu_schreier_vec[0] = 0;
	      orbit[i]->trans = NULL;
	   }
				/*schreier_vec[0] = Laenge des Wortes*/
  erz_anz = Stab->gen_no;
  Ssize = 0;
  Size = MIN_SPEICHER;
  histsize = 0;
  histinvsize = 0;
  histanz = 0;
  histinvanz = 0;
  *leng = 1;
  old = 0;
  count = *leng;
  Sverz = NULL;
  Erz = (matrix_TYP **) malloc(erz_anz *sizeof(matrix_TYP *));
  for(i=0;i< Stab->gen_no;i++)
    Erz[i] = Stab->gen[i];
  
  Erz_inv = (matrix_TYP **) malloc(erz_anz *sizeof(matrix_TYP *));
  for(i=0;i<erz_anz;i++)
       Erz_inv[i] = mat_inv(Erz[i]);

  I = einheitsmatrix(Stab->gen[0]->rows);
  I1 = einheitsmatrix(Stab->gen[0]->rows);
  I2 = einheitsmatrix(Stab->gen[0]->rows);
  matrix_speichern(I1,&hist,&histsize,&histanz);
  matrix_speichern(I2,&histinv,&histinvsize,&histinvanz);
  orbit[0]->elm = copy_mat(Y); 

  base = get_base(Stab);
  strong = NULL;

      merk_stab = (matrix_TYP **)malloc(sizeof(matrix_TYP*)); 

  re.z = -1; re.n = 1; li.z = 1; li.n = 1;
  neu_anz = 1;
  while(neu_anz !=0){
     neu_anz = 0;

     for(i=old; i<*leng; i++){
        for(j=0; j<erz_anz; j++){ 
           test = mat_mul(Erz[j],orbit[i]->elm);
           Check_mat(test);
           schonda = is_neu(test,orbit,count); 
           if(schonda == -1){
             Size = append_new_element(&orbit,count,test,i,j,
                                Stab,SW,Mat,Mat_inv, matanz,Size);
          
             K1 = mat_mul(Erz[j],hist[i]);
             matrix_speichern(K1, &hist, &histsize, &histanz);
             K1 = NULL;
             K1 = mat_mul(histinv[i], Erz_inv[j]);
             matrix_speichern(K1, &histinv, &histinvsize, &histinvanz);
             K1 = NULL;
         
              count++;
              neu_anz++;
           }/*schonda*/
           else{ 
		   /*if(schonda != -1) schonda = Nr des gleichen El.*/
             		
             K1 = mat_mul(Erz[j], hist[i]);
             K2 = mat_mul(histinv[schonda], K1);

	     if(mat_vergleich(K2,I) != 0){  
	      Liste(K2, &neu_Stab->gen, &neu_Stab->gen_no, &sch, &Ssize);
              if(sch == -1){
 	      old_gen_no = neu_Stab->gen_no;
 	      red = red_gen(neu_Stab,base,&strong,neu_Stab->gen_no-1);

		/* neues Stabilisatorelement wurde genommen: */

	         if(neu_Stab->gen_no == old_gen_no){

                   /*erzeuge Wort, welches das Stabilisator-Element in den 
                     Erz der Gr. (also des vorhergehenden Stabilisators) 
                     ausdrueckt*/

	           Word(j,i,neu_SW, neu_Stab->gen_no, orbit, schonda);

	           merk_stab[neu_Stab->gen_no -1] =  copy_mat(K2);
                   merk_stab = (matrix_TYP **)realloc(merk_stab,
                             (neu_Stab->gen_no+1)*sizeof(matrix_TYP *));

                 }/* if(neu_Stab->gen_no != old_gen_no) **/

                /*!Permutation der Stabilisator-Elemente 
                noch aendern (rueckgaengig machen)!! */

		for(l=0; l<neu_Stab->gen_no; l++){
		   if(mat_comp(merk_stab[l],neu_Stab->gen[l]) != 0){
			ll = Finde(merk_stab[l], neu_Stab,l);
			tmp = copy_mat(neu_Stab->gen[ll]);
			free_mat(neu_Stab->gen[ll]); neu_Stab->gen[ll] = NULL;
			neu_Stab->gen[ll] = copy_mat(neu_Stab->gen[l]);
			free_mat(neu_Stab->gen[l]); neu_Stab->gen[l] = NULL;
			neu_Stab->gen[l] = copy_mat(tmp);
			free_mat(tmp);	tmp = NULL;
		   }
     		}

	      }/**if(sch == -1) **/
             }/** if( K2 != I) **/
           if(schonda == -1 || (schonda != -1 && sch != -1)){
              free_mat(K2); K2 = NULL;
           }
           free_mat(K1); K1 = NULL;  K2 = NULL;

           free_mat(test); test = NULL;
	   }/**else **/
        }/*for(j)*/
     }/*for(i)*/

     old = *leng;
     *leng = count;
  }/*while*/

    for(i=0; i<histanz; i++)
    {
      free_mat(hist[i]); free_mat(histinv[i]);
    }
    free(hist); free(histinv);
    free_mat(I); I = NULL;

  for(i=0;i<erz_anz;i++){
      free_mat(Erz_inv[i]); Erz_inv[i] = NULL;
  }
  free(Erz); free(Erz_inv);
  
  for(i=(*leng); i<Size; i++){
     /* free(orbit[i]->schreier_vec);
     free(orbit[i]->neu_schreier_vec); */
     free(orbit[i]);
  }

orbit=(neu_ele_TYP  **)realloc(orbit,(*leng) *sizeof(neu_ele_TYP *));
erg = (ele_TYP **)malloc((*leng)*sizeof(ele_TYP*));
for(i=0; i<(*leng); i++){
   erg[i]= (ele_TYP*)malloc(sizeof(ele_TYP));
   erg[i]->elm = copy_mat(orbit[i]->elm);
   free_mat(orbit[i]->elm); orbit[i]->elm = NULL;
   erg[i]->schreier_vec = 
              (int*)malloc((orbit[i]->neu_schreier_vec[0]+1)*sizeof(int));
   memcpy(erg[i]->schreier_vec, orbit[i]->neu_schreier_vec, 
                           (orbit[i]->neu_schreier_vec[0]+1)*sizeof(int));
   free(orbit[i]->schreier_vec);
   free(orbit[i]->neu_schreier_vec);
   erg[i]->trans = copy_mat(orbit[i]->trans);
   if(orbit[i]->trans != NULL){
      free_mat(orbit[i]->trans); orbit[i]->trans = NULL;
   }
}
free(orbit);

/*umschreiben von neu_SW auf Erzeuger der Raumgruppe*/

change(neu_SW, SW, neu_Stab->gen_no,Stab, Mat, Mat_inv,matanz);

return(erg);
}/* stab_orbit() */
