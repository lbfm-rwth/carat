
/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: barzen.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/

#include "typedef.h"
#include "matrix.h"
#include "getput.h"
#include "tietzetrans.h"
#include "presentation.h"
#include "sort.h"

static derivedsg_TYP * make_dersg(int gen_no);

/* convert word_TYP to int* */
/**************************************************************************\
@-------------------------------------------------------------------------
@int *make_prod(word_TYP *word)
@
@converts a word_TYP to a list of integers representing a word in the given
@generators
@
@-------------------------------------------------------------------------
\**************************************************************************/
int *make_prod(word_TYP *word)
{
 int  	i,j;
 int	*erg;
 int	a,b,c;

 a = 0;
 for(i=0; i<word->trans->rows-1; i++){
    if(word->trans->kgv != 1) exit(3);
    a += abs(word->trans->array.SZ[i][0]);
 }
 b = a + word->word[0];
 erg = (int*)calloc(b+1,sizeof(int));

 c = 1;
 for(j=0; j<word->trans->rows-1; j++){
    if(word->trans->array.SZ[j][0] > 0){
    for(i=0; i< abs(word->trans->array.SZ[j][0]); i++)
    erg[c++] = j+1;
    }/*if ()*/
    else{
    for(i=0; i< abs(word->trans->array.SZ[j][0]); i++)
    erg[c++] = (-1)*(j+1);
    }/*else()*/
 }/* for() */

 for(j=0; j<word->word[0]; j++)
    erg[c++] = word->word[j+1];

/* for(j=a;j<b; j++)
    erg[c++] = word->word[0];*/

 erg[0] = c-1;

 return(erg);
}/** make_prod() **/

/**************************************************************************\
@--------------------------------------------------------------------------
@static derivedsg_TYP * make_dersg(int gen_no)
@
@allocates memory for a derivedsg_TYP 
@and sets ERG->firstfree = gen_no.
@
@--------------------------------------------------------------------------
\**************************************************************************/
static derivedsg_TYP * make_dersg(int gen_no)
{
 int		i;
 derivedsg_TYP *ERG;

	ERG = (derivedsg_TYP *)malloc(sizeof(derivedsg_TYP));
	ERG->list = (derived_TYP **)malloc((EXT_SIZE) *sizeof(derived_TYP *));
	ERG->firstfree = gen_no;
	ERG->sizemult = 1;
	for(i=0; i<gen_no; i++){
	   ERG->list[i] = (derived_TYP *)malloc(sizeof(derived_TYP));
	   ERG->list[i]->element = NULL; 
	   ERG->list[i]->product = NULL; 
	   ERG->list[i]->nproduct =0;
	   ERG->list[i]->left = NIL;
	   ERG->list[i]->right = NIL;
 	}
 return(ERG);
}/** make_dersg() **/

 
/**************************************************************************\
@--------------------------------------------------------------------------
@ presentation_TYP* make_pres(anne_presentation_TYP *PR)
@
@transfers anne_presentation_TYP *PR to presentation_TYP, *PR is set free 
@--------------------------------------------------------------------------
\**************************************************************************/
presentation_TYP* make_pres(anne_presentation_TYP *PR)
{
 presentation_TYP *erg;
 int i;
 
 erg = (presentation_TYP *)malloc(sizeof(presentation_TYP));
 erg->generators = make_dersg(PR->gen_no); 
 erg->relators = (relator_TYP *)malloc((EXT_SIZE) *sizeof(relator_TYP));
 erg->norelators = PR->rel_no;
 /*erg->norelators = PR->gen_no + PR->rel_no; */
 erg->ext_factor = erg->norelators;

for(i=0; i<PR->gen_no; i++){
   erg->generators->list[i]->element = PR->generators[i];
   erg->generators->list[i]->nproduct = 0;
   erg->generators->list[i]->product = NULL;
   erg->generators->list[i]->left = 0;
   erg->generators->list[i]->right = 0;
}

for(i=0; i<PR->rel_no; i++){
   erg->relators[i].lhsproduct = NULL;
   erg->relators[i].lhsnproduct = 0; 
   erg->relators[i].rhsproduct = PR->relators[i].lhsproduct;
   erg->relators[i].rhsnproduct = PR->relators[i].lhs_no;
}
 free(PR);
return(erg);
}/** make_pres() **/

/**************************************************************************\
@--------------------------------------------------------------------------
@matrix_TYP *barzen(polyeder_TYP *Pol)
@
@returns a matrix which is the baryzentrum of Pol
@
@--------------------------------------------------------------------------
\**************************************************************************/
matrix_TYP *barzen(polyeder_TYP *Pol)
{
  matrix_TYP *erg;
  int	     *list;
  int	     i,j,d,kgv;

  d = Pol->vert[0]->dim;
  erg = init_mat(d,1,"");

	/*bestimme Baryzentrum der Wand*/
      list = (int*) calloc(Pol->vert_no,sizeof(int));
      for(i=0; i<Pol->vert_no; i++) 		
      {
      list[i] = Pol->vert[i]->v[d-1];
      }
      kgv = KKGV(list, Pol->vert_no);	/*funktion KKGV extern */ 
  
      for(j=0; j<d; j++){
          for(i=0; i<Pol->vert_no; i++)
          {
          	erg->array.SZ[j][0] +=  (kgv/list[i])*(Pol->vert[i]->v[j]); 
          }
      }
      erg->kgv = erg->array.SZ[d-1][0];
      Check_mat(erg);
    free(list);
return(erg);
}/** barzen() **/

/*
    Liste[i] = j -> P->wall[i] hat Seitentrafo = P->wall[j-1]->mat falls j>0 
			     und Seitentrafo = (P->wall[j-1]->mat)^-1 ; j<0.
*/
/**************************************************************************\
@--------------------------------------------------------------------------
@ int* generate_Liste(polyeder_TYP *Pol,int **Liste_inv)
@
@ creates a list of integers decoding the sidepairing of a given polyeder Pol
@ with sidetransformations (in Pol->mat).  
@ List_inv will contain information of the inverse sidetrasformations
@--------------------------------------------------------------------------
\**************************************************************************/
int* generate_Liste(polyeder_TYP *Pol,int **Liste_inv)
{
 int 			*Liste;
 int 			i,j,s,a;
 matrix_TYP 		*M_inv;
 int 			*test;

    test = (int*)malloc(Pol->wall_no * sizeof(int));
    for(i=0; i<Pol->wall_no; i++)
    test[i] = 1;

    Liste = (int*)malloc(Pol->wall_no * sizeof(int));
    *Liste_inv = (int*)malloc(Pol->wall_no * sizeof(int));

    for(i=0; i<Pol->wall_no; i++){
	if(test[i] != 0){
         Check_mat(Pol->wall[i]->mat);
	 M_inv = mat_inv(Pol->wall[i]->mat);
	    a = 0;
            for(j=i; j<Pol->wall_no; j++){
          	if(test[j]!=0){
	        s = mat_comp(Pol->wall[j]->mat,M_inv);
	          if(s == 0){
		   test[i]=0; test[j]=0;
		   Liste[j] = (-1)*(i+1);
		   Liste[i] = (i+1);
		   Liste_inv[0][i] = j;
		   Liste_inv[0][j] = i;
		   a = 1;
                   j = Pol->wall_no;
                  }
               }
            }/*for(j)*/
	    if(a == 0){
	    printf("wrong sidetransformations given \n");
	    exit(3);
	    }
	 free_mat(M_inv); M_inv = NULL;
        }
    }/*for(i)*/
    free(test);

return(Liste);
}/** generate_Liste(Pol) **/

