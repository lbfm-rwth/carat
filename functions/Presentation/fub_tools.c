
/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: fub_tools.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

#include "typedef.h"
#include "matrix.h"
#include "symm.h"
#include "getput.h"
#include "bravais.h"
#include "polyeder.h"
#include "presentation.h"
#include "tools.h"

/**************************************************************************\
@--------------------------------------------------------------------------
@ matrix_TYP *word_to_mat(int *word, matrix_TYP **Mat, 
@ 			matrix_TYP **Mat_inv, int matanz)
@
@ calculates the matrix determined by word (multiplying the matrices in Mat 
@ and Mat_inv according to the numbers occuring in the integer string word)
@
@--------------------------------------------------------------------------
\**************************************************************************/
matrix_TYP *word_to_mat(int *word, matrix_TYP **Mat, 
			matrix_TYP **Mat_inv, int matanz)
{
 matrix_TYP 	*A, *B;
 int		i;
 
 A = init_mat(Mat[0]->rows, Mat[0]->cols, "1");

 for(i=1; i<= word[0]; i++){
    if(word[i]>0)
      B = mat_mul(A,Mat[(word[i])-1]);
    else
      B = mat_mul(A,Mat_inv[-(word[i])-1]);
    free_mat(A); A = NULL;
    A = copy_mat(B);
    free_mat(B); B = NULL;
 }/**for(i)**/

 return(A);
}/* function word_to_mat */

 /* noch das Element Orbit[i]->trans als Translation an word_to_mat anfuegen */
/**************************************************************************\
@--------------------------------------------------------------------------
@ matrix_TYP *orb_word_to_mat(int *word, matrix_TYP **Mat, 
@ 			matrix_TYP **Mat_inv, int matanz, matrix_TYP *trans)
@
@ calculates the matrix given by a word and a translation (in case of affine
@ groups), Mat and Mat_inv contain the generators of a group and their inverses 
@
@--------------------------------------------------------------------------
\**************************************************************************/
 matrix_TYP *orb_word_to_mat(int *word, matrix_TYP **Mat, 
 			matrix_TYP **Mat_inv, int matanz, matrix_TYP *trans)
 {
   int  		i;
   matrix_TYP 	*A;
 
    A = word_to_mat(word, Mat, Mat_inv, matanz);
    for(i=0; i<A->rows-1; i++) 
    A->array.SZ[i][A->cols-1] += (A->kgv)*trans->array.SZ[i][0];
 
  return(A);
 }/** orb_word_to_mat() **/
 
/**************************************************************************\
@--------------------------------------------------------------------------
@ matrix_TYP *Stab_word_to_mat(Stab_word_TYP SW, matrix_TYP **Mat, 
@ 			     matrix_TYP **Mat_inv, int matanz, matrix_TYP *Tmat)
@
@ transforms a Stab_word_TYP SW into a matrix 
@ Mat and Mat_inv contain the generators of a space group and their inverses 
@ the rows of Tmat determine the translation lattice.
@--------------------------------------------------------------------------
\**************************************************************************/
matrix_TYP *Stab_word_to_mat(Stab_word_TYP SW, matrix_TYP **Mat, 
			     matrix_TYP **Mat_inv, int matanz, matrix_TYP *Tmat)
{
 matrix_TYP *A, *B, *C;
 int i,s;

 s = Tmat->rows; 
 A = word_to_mat(SW.sword, Mat, Mat_inv, matanz);

 B = copy_mat(SW.trans);

 C = init_mat(s+1, s+1, "1");
     for(i=0; i< SW.trans->rows; i++){
         C->array.SZ[i][i] = B->kgv;
	 C->array.SZ[i][s] = B->array.SZ[i][0];
        }

	 free_mat(B); B=NULL;
 B = mat_mul(C,A);
 	free_mat(A); A=NULL;
 	free_mat(C); C=NULL;

return(B);
}/***Stab_word_to_mat()***/

/* bilde Ungleichung bAb-mAm+2*(m-b)A > 0 , also berechne bAb-mAm-2*(b-m)A 
   wobei m der Startvektor und b = gm+t.
 */
/**************************************************************************\
@--------------------------------------------------------------------------
@ int *wand_ungl(matrix_TYP *M, matrix_TYP *B,matrix_TYP *Form)
@
@ calculates the wall-inequality between the points  M and B according to the 
@ metric given by Form. 
@
@--------------------------------------------------------------------------
\**************************************************************************/
int *wand_ungl(matrix_TYP *M, matrix_TYP *B,matrix_TYP *Form)
{
    int 	i, g;
    matrix_TYP  *k1, *k2;
    rational 	li, re, kk;
    matrix_TYP 	*AA,*A;
    int 	*erg;

    A = tr_pose(M);
    A->cols --;
    k1 = scal_pr(A, Form, TRUE);
    if(k1->kgv == 0)
       k1->kgv = 1;
    free_mat(A); A = NULL;
    A = tr_pose(B);
    A->cols --;
    k2 = scal_pr(A, Form, TRUE);
    if(k2->kgv == 0)
       k2->kgv = 1;
    free_mat(A); A = NULL;

    if(k1->kgv != k2->kgv){
        kk.z = (k1->kgv)*(k2->array.SZ[0][0]) - (k2->kgv)*(k1->array.SZ[0][0]);
        kk.n = (k1->kgv)*(k2->kgv);
        g = GGT(kk.z,kk.n);
        kk.z = kk.z/g;
        kk.n = kk.n/g;
    }
    else{
        kk.z = k2->array.SZ[0][0] - k1->array.SZ[0][0];
        kk.n = k1->kgv;
        g = GGT(kk.z,kk.n);
        kk.z = kk.z/g;
        kk.n = kk.n/g;
    }

    li.z = -1;
    li.n = 1;
    re.z = 1;
    re.n = 1;
    AA = mat_add(B, M, li, re);   /* M-B */   /*anne, 14.8.*/
    AA->rows --;
    A = mat_mul(Form, AA);
    AA->rows ++;
    free_mat(AA); AA = NULL;
    sscal_mul(A,2);		/* anne, 14.8.*/
    Check_mat(A);
    erg = (int *) malloc(M->rows * sizeof(int));
    if(kk.n > 0){
       for (i = 0; i < M->rows - 1; i++)
   	   erg[i] = (kk.n)*A->array.SZ[i][0];
    erg[M->rows - 1] = kk.z * A->kgv;		/* anne, 15.8. *A->kgv */
    }
    else{
       for (i = 0; i < M->rows - 1; i++)
   	   erg[i] = (-1)*(kk.n)*A->array.SZ[i][0];
    erg[M->rows - 1] = (-1)*kk.z;
    }
    free_mat(k1); k1 = NULL;
    free_mat(k2); k2 = NULL;
    free_mat(A); A = NULL;
return(erg);
}/***wand_ungl()***/

/*erzeugt Waende zwischen orbit[0] und allen weiteren Elementen orbit[i],
  schneidet die relevanten Waende vom Polyeder ab und gibt alle Waende wieder
  frei */
/**************************************************************************\
@--------------------------------------------------------------------------
@ int wand()
@
@ int erz_anz                  number of goup generators
@ matrix_TYP **Erz, **Erz_inv  generators of a group and their inverses
@ ele_TYP ** orbit             orbit of a point orbit[0] under the group 
@                              generated by Erz
@ int length                   length of the orbit
@ polyeder_TYP * Pol           a given polyeder, orbit[0] inner point of Pol
@ matrix_TYP * Form            positive definite quadratic form according to 
@                              which Pol is given
@ int **rel_kand               a string of 0 and 1, rel_kand[0] = 1 and
@                              rel_kand[i]=1 if the wall between orbit[i] and
@                              orbit[0] was used to refine Pol, 0 else.     
@
@ calculates the walls between orbit[0] and every other orbit[i]
@ and cuts of these walls from the polyeder Pol. All calculated walls
@ are set free afterwards. 
@ The function returns the number of refinements made on Pol.
@
@--------------------------------------------------------------------------
\**************************************************************************/
int wand(int erz_anz, matrix_TYP **Erz, matrix_TYP **Erz_inv, ele_TYP ** orbit, int length, polyeder_TYP * Pol, matrix_TYP * Form, int **rel_kand)
{
    int i, j;
    int weiter;
    int dim;
    wall_TYP **wall;
    int *erg;

    dim = orbit[0]->elm->rows;
    erg = (int*) calloc(length,sizeof(int));
    wall = (wall_TYP **) malloc((length - 1) * sizeof(wall_TYP *));
    for (i = 0; i < length - 1; i++) {
	wall[i] = init_wall(dim);
	wall[i]->mat = orb_word_to_mat(orbit[i + 1]->schreier_vec,
			     Erz, Erz_inv, erz_anz, orbit[i + 1]->trans);
        Check_mat(wall[i]->mat);
        if(wall[i]->gl != NULL) free(wall[i]->gl);
	wall[i]->gl = wand_ungl(orbit[0]->elm, orbit[i + 1]->elm, Form);
   	wall[i]->word = (word_TYP *)malloc(sizeof(word_TYP));
 	wall[i]->word->dim = dim;
        wall[i]->neu = 0;
        wall[i]->next_no = 0;
        wall[i]->next = NULL;
        wall[i]->ext_no = 0;
        wall[i]->extra = NULL;
	wall[i]->word->trans = copy_mat(orbit[i+1]->trans);
        wall[i]->word->word = 
		(int*)malloc((orbit[i + 1]->schreier_vec[0]+1) * sizeof(int));
	memcpy(wall[i]->word->word,orbit[i + 1]->schreier_vec,
			(orbit[i + 1]->schreier_vec[0]+1) * sizeof(int));
/* anne, 17.8.
        wall[i]->word->word = 
		(int*)malloc((orbit[i + 1]->schreier_vec[0]) * sizeof(int));
	memcpy(wall[i]->word->word,orbit[i + 1]->schreier_vec,
			(orbit[i + 1]->schreier_vec[0]) * sizeof(int));
*/
    }
    weiter = 0;
    for (i = 0; i < length - 1; i++) {
        normal_wall(wall[i]);
	j = refine_polyeder(Pol, wall[i]);
        if(j == 1){ 		/* wall[i] ist relevante Wand */
           erg[i+1] = 1;	/* wall[i] = wand zwischen orb[0] & orb[i+1] */
           weiter ++;
        }
	free_wall(&wall[i]);
    }
    erg[0] = 1;
    *rel_kand = erg;
    free(wall);
return(weiter);
}/**wand()**/

/**************************************************************************\
@--------------------------------------------------------------------------
@ void P_ausgabe(polyeder_TYP *Pol)
@
@ prints the coordinates of the polyeder Pol and its sidetransformations 
@ and the words how the sidetransformations are expressed in the spacegroup.
@ 
@--------------------------------------------------------------------------
\**************************************************************************/
void P_ausgabe(polyeder_TYP *Pol)
{
  int	i,j;
put_polyeder(Pol);
for(i=0; i<Pol->wall_no; i++){
   put_mat(Pol->wall[i]->mat,NULL,"",0);
   put_mat(Pol->wall[i]->word->trans,NULL,"word_trans",0);
   for(j=0; j<=Pol->wall[i]->word->word[0]; j++)
       printf("%d " ,Pol->wall[i]->word->word[j]);
   printf("\n");
  }
}/** P_ausgabe() **/

/* Funktion welche zwei Stabilisator_Woerter miteinander ausmultipliziert

merke mir worte als (vektor);wort(Nr. der Erzeuger)
wobei vector die Koordinaten eines Translationsvektors
enthaelt.
Aufruf:
Stab_word_TYP erg;
Stab_word_mul(sw1, sw2, Mat, Mat_inv, matanz, &erg);

typedef struct {
int **sword;
matrix_TYP *trans;
}Stab_word;

schreibe funktion : multipliziere zwei worte
(vec1);(wort1)*(vec2);(wort2) = (vec1+vec3);(wort1)*(wort2)
wobei vec3 = wort1 auf vec2 angewandt ist. 
(geht am leichtesten indem man wort1 als matrix ausrechnet, 
 und den LINEAREN ANTEIL dann mit vec2 multipliziert)
wort1*wort2 ist einfach das aneinanderhaengen der worte.
 */

/**************************************************************************\
@--------------------------------------------------------------------------
@ void Stab_word_mul(sw1, sw2, Mat, Mat_inv, matanz, erg)
@ Stab_word_TYP *sw1, *sw2;
@ matrix_TYP **Mat, **Mat_inv; 
@ int matanz;
@ Stab_word_TYP *erg;
@
@ multiplies two words (Stab_word_TYP) sw1 and sw2 and puts the result to erg.
@ Mat and Mat_inv are the matrices (group generators and inverses) according
@ to which the words are build. matanz is the number of matrices Mat.
@
@--------------------------------------------------------------------------
\**************************************************************************/
void Stab_word_mul(sw1, sw2, Mat, Mat_inv, matanz, erg)
Stab_word_TYP *sw1, *sw2;
matrix_TYP **Mat, **Mat_inv; 
int matanz;
Stab_word_TYP *erg;
{
 int 		i;
 int 		a,c;
 rational	re,li;
 matrix_TYP     *A,*B;

 a = (*sw1).sword[0]+(*sw2).sword[0];
 (*erg).sword = (int *)malloc((a+1)*sizeof(int));
 (*erg).sword[0] = a;
 for(i=0,c=1; i<(*sw1).sword[0]; i++){
     (*erg).sword[c] = (*sw1).sword[c];
     c++;
    }
 for(i=0; i<(*sw2).sword[0]; i++){
     (*erg).sword[c++] = (*sw2).sword[i+1];
    }
 
 A = word_to_mat((*sw1).sword, Mat, Mat_inv, matanz);
/* wende nur den linearen Anteil von A auf (*sw2).trans an, setze letzte Spalte von A = 0 !*/

 for(i=0; i<A->rows-1; i++)
	A->array.SZ[i][A->cols-1] = 0;

 if((*sw2).trans != NULL)
    B = mat_mul(A,(*sw2).trans);
 else{
 /* (*sw2).trans ist der Null-vektor */
    (*sw2).trans = init_mat(A->rows, 1,"");
    (*sw2).trans->array.SZ[A->rows-1][0] = 1;

    B = mat_mul(A,(*sw2).trans);
 }
/*
put_mat(A,NULL, "A", 0);
put_mat(B,NULL, "B", 0);
*/
free_mat(A); A = NULL;

li.z = B->kgv; li.n = 1; re.z = (*sw1).trans->kgv; re.n = 1;

if((*sw1).trans != NULL)
A = mat_add((*sw1).trans,B,li,re);
else{
    (*sw1).trans = init_mat(B->rows, 1,"0");
    (*sw1).trans->array.SZ[B->rows-1][0] = 1;
    A = mat_add((*sw1).trans,B,li,re);
}
A->array.SZ[A->rows-1][0] = (B->kgv)*((*sw1).trans->kgv);
Check_mat(A);
/*put_mat(A,NULL, "A", 0);*/
 
 (*erg).trans = A;
free_mat(B); B = NULL;
/*
 for(i=0; i<a; i++){
     printf("%d ",(*erg).sword[i]); 
     printf("\n");
    }
*/
}/** Stab_word_mul() **/

/**************************************************************************\
@---------------------------------------------------------------------------
@ void free_ele(orb)
@ ele_TYP *orb; 
@ 
@    frees, what was allocated in *orb and sets (*orb) = NULL 
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
void free_ele(orb)
ele_TYP *orb; 
{
   if (orb != NULL){
      if (orb->elm != NULL){
         free_mat(orb->elm); orb->elm = NULL;
      }
      if (orb->trans != NULL){
         free_mat(orb->trans); orb->trans = NULL;
      }
/*      if (orb->schreier_vec[0] > 0 && orb->schreier_vec != NULL)  */
      if (orb->schreier_vec != NULL)
         free(orb->schreier_vec);
      free(orb);
      orb = NULL;
   }
}/** free_ele() **/

/**************************************************************************\
@---------------------------------------------------------------------------
@ void free_Stab_word(SW)
@ Stab_word_TYP *SW;
@ 
@    frees, what was allocated in *SW and sets SW = NULL;
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
 void free_Stab_word(SW)
 Stab_word_TYP *SW;
{
  if((*SW).sword[0] != 0 && (*SW).sword[0] != NULL)
     free((*SW).sword);
  if((*SW).trans != NULL){
     free_mat((*SW).trans); (*SW).trans = NULL;
   }
  free(SW);
  (SW) = NULL;
}/** free_Stab_word() **/

/**************************************************************************\
@---------------------------------------------------------------------------
@ void put_Stab_word(SW)
@ Stab_word_TYP *SW;
@ 
@    puts, what was allocated in *SW 
@
@---------------------------------------------------------------------------
\**************************************************************************/
 void put_Stab_word(SW)
 Stab_word_TYP *SW;
{
  int i;
  if((*SW).sword[0] != 0 && (*SW).sword[0] != NULL)
  for(i=0; i<=(*SW).sword[0]; i++)
      printf("%d ", (*SW).sword[i]);
  printf("\n");
  if((*SW).trans != NULL)
     put_mat((*SW).trans,NULL,"SW",0);
}/* put_Stab_word */

/**************************************************************************\
@--------------------------------------------------------------------------
@ void Word_mul(sw1, sw2, Mat, Mat_inv, matanz, erg)
@ word_TYP *sw1, *sw2;
@ matrix_TYP **Mat, **Mat_inv; 
@ int matanz;
@ word_TYP *erg;
@
@ multiplies two words sw1 and sw2 and puts the result to erg.
@ Mat and Mat_inv are the matrices (group generators and inverses) according
@ to which the words are build. matanz is the number of matrices Mat.
@--------------------------------------------------------------------------
\**************************************************************************/
void Word_mul(sw1, sw2, Mat, Mat_inv, matanz, erg)
word_TYP *sw1, *sw2;
matrix_TYP **Mat, **Mat_inv; 
int matanz;
word_TYP *erg;
{
 int 		i;
 int 		a,c;
 rational	re,li;
 matrix_TYP     *A,*B;

 a = (*sw1).word[0]+(*sw2).word[0];
 (*erg).word = (int *)malloc((a+1)*sizeof(int));
 (*erg).word[0] = a;
 for(i=0,c=1; i<(*sw1).word[0]; i++){
     (*erg).word[c] = (*sw1).word[c];
     c++;
    }
 for(i=0; i<(*sw2).word[0]; i++){
     (*erg).word[c++] = (*sw2).word[i+1];
    }
 
 A = word_to_mat((*sw1).word, Mat, Mat_inv, matanz);
 if((*sw2).trans != NULL)
    B = mat_mul(A,(*sw2).trans);
 else{
 /* (*sw2).trans ist der Null-vektor */
    (*sw2).trans = init_mat(A->rows, 1,"");
    (*sw2).trans->array.SZ[A->rows-1][0] = 1;

    B = mat_mul(A,(*sw2).trans);
 }
/*
put_mat(A,NULL, "A", 0);
put_mat(B,NULL, "B", 0);
*/
free_mat(A); A = NULL;

li.z = B->kgv; li.n = 1; re.z = (*sw1).trans->kgv; re.n = 1;

if((*sw1).trans != NULL)
A = mat_add((*sw1).trans,B,li,re);
else{
    (*sw1).trans = init_mat(B->rows, 1,"0");
    (*sw1).trans->array.SZ[B->rows-1][0] = 1;
    A = mat_add((*sw1).trans,B,li,re);
}
A->array.SZ[A->rows-1][0] = (B->kgv)*((*sw1).trans->kgv);
Check_mat(A);
/*put_mat(A,NULL, "A", 0);*/
 
 (*erg).trans = A;
free_mat(B); B = NULL;
/*
 for(i=0; i<a; i++){
     printf("%d ",(*erg).word[i]); 
     printf("\n");
    }
*/
}/** Word_mul() **/
