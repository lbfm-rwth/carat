#include"typedef.h"
#include"getput.h"
#include"bravais.h"
#include"matrix.h"
#include"tietzetrans.h"

int INFO_LEVEL;
int SFLAG;

/****************************************************************************
@
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@ FILE: tietze.c
@----------------------------------------------------------------------------
@----------------------------------------------------------------------------
@
*****************************************************************************/



/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ void free_presentation(presentation_TYP *p)
@
@ frees all records of p, and p itself.
@ p will also be changed.
@
@----------------------------------------------------------------------------
@
*****************************************************************************/
void free_presentation(presentation_TYP *p)
{

    int i;

    relator_TYP *rel=p->relators;

    derivedsg_TYP *gen=p->generators;

    for (i=0;i<p->norelators;i++){
       if (rel[i].lhsnproduct > 0 && rel[i].lhsproduct != NULL){
          free(rel[i].lhsproduct);
       }
       if (rel[i].rhsnproduct > 0 && rel[i].rhsproduct != NULL){
          free(rel[i].rhsproduct);
       }
    }

    if (p->norelators >0  && rel != NULL){
       free(rel);
       p->relators = NULL;
    }

    p->generators = NULL;
    if (gen != NULL){

       for (i=0;i<gen->firstfree;i++){
          if (gen->list[i] != NULL){
             if (gen->list[i]->element != NULL){
                free_mat(gen->list[i]->element);
             }
             if (gen->list[i]->nproduct > 0 && gen->list[i]->product != NULL){
                free(gen->list[i]->product);
             }
             free(gen->list[i]);
          }
       }

       if (gen->list != NULL){
          free(gen->list);
       }
       free(gen);
    }

    free(p);

    return;
}

/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ int simplify_relator(relator_TYP *r)
@
@ simplify_relator conjugates the relator r in the free group to make
@ it shorter. It  might invert the relator if this leads to less ^(-1)
@ in total. The return is TRUE iff the resulting relator has length 0.
@ The function does not "cyclicly reduce" the relators.
@
@ CAUTION: the function assumes the relator to be totaly on the rigth
@ hand side of the equation.
@ NOTE: The function will not realloc/free any of the components of r.
@
@----------------------------------------------------------------------------
@
*****************************************************************************/
int simplify_relator(relator_TYP *r)
{
   int i,
       j,
       null_flag = (r->rhsnproduct == 0),
       recurse_flag = FALSE;

   /* check trivialities */
   if (r->lhsnproduct != 0){
      fprintf(stderr,"simplify_relator should only be used for relators\n");
      fprintf(stderr,"without left hand side\n");
      exit(3);
   }

   /* firstly scan the relator for occurences of xi*xi^-1 or xi^-1*xi */
   for (i=1;i<r->rhsnproduct;i++){
      if (r->rhsproduct[i-1] == (-r->rhsproduct[i])){
         recurse_flag = TRUE;
         r->rhsnproduct -= 2;
         for (j=i-1;j<r->rhsnproduct;j++){
            r->rhsproduct[j] = r->rhsproduct[j+2];
         }
      }
   }

   /* secondly check whether this relator is conjugate to a shorther one
      in the free group */
   if (r->rhsnproduct > 1 &&
      (r->rhsproduct[0] == (-r->rhsproduct[r->rhsnproduct-1]))){
      r->rhsnproduct -= 2;
      for (i=0;i<r->rhsnproduct;i++){
         r->rhsproduct[i] = r->rhsproduct[i+1];
      }
      recurse_flag = TRUE;
   }

   if (recurse_flag && r->rhsnproduct > 0){
      simplify_relator(r);
   }

   /* now basicaly everything is done, but it we migth be better off with
      the inverse of this generator */
   j = 0;
   for (i=0;i<r->rhsnproduct;i++){
      if (r->rhsproduct[i] > 0){
         j++;
      }
      else if (r->rhsproduct[i] < 0){
         j--;
      }
      else{
         fprintf(stderr,"error in simplify_relator\n");
         exit(3);
      }
   }

   /* if there are more ^(-1) then positive signs, invert the relator */
   if (j<0){
      for (i=0;2*i<r->rhsnproduct;i++){
         j = -r->rhsproduct[r->rhsnproduct-i-1];
         r->rhsproduct[r->rhsnproduct-i-1] = -r->rhsproduct[i];
         r->rhsproduct[i]=j;
      }
   }

   if (r->rhsnproduct == 0 && !null_flag && r->rhsproduct != NULL){
      free(r->rhsproduct);
      r->rhsproduct = NULL;
   }

   return (r->rhsnproduct == 0);
}

/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ void put_presentation(presentation_TYP *pres,char *file,char *option)
@
@ puts the presentation in pres. It can handle gap format and the carat
@ format. It also chechs the presentation for illegal "0" in one of the
@ given relators.
@
@ presentation_TYP *pres: The presentation in question. All relators have
@                         to have an empty left hand side, otherwise
@                         the function will exit with a message to stderr.
@ char *file            : The name of the output file. If file == NULL all
@                         output will be redirected to stdout.
@ char *option          : Drives the behaviour of the function. If
@                         option == "G", gap code is produced,
@                         option == "C", carat code is produced.
@ 
@
@----------------------------------------------------------------------------
@
*****************************************************************************/
void put_presentation(presentation_TYP *pres,char *file,char *option)
{
   int i,
       j,
       max=0;

   relator_TYP *rel=pres->relators;

   FILE *OUT;

   /* open the outfile if nessesary */
   if (file == NULL){
      OUT = stdout;
   }
   else{
      OUT = fopen(file,"w");
   }

   if (strcmp("G",option) == 0){

     /* The GAP part */

     fprintf(OUT,"F:=FreeGroup(%d);\n",pres->generators->firstfree);

     /* output definitions x1:=F.1; .... */
     for (i=0;i<pres->generators->firstfree;i++){
         fprintf(OUT,"x%d := F.%d;\n",i+1,i+1);
     }

     fprintf(OUT,"H:=F/[\n");
     for (i=0;i<pres->norelators;i++){

        /* don't handle the case of a left hand side jet */
        if (rel[i].lhsnproduct > 0){
           fprintf(stderr,"error in put_presentation\n");
           exit(3);
        }

        for (j=0;j<rel[i].rhsnproduct;j++){
           if (j != 0) fprintf(OUT,"*");
           if (rel[i].rhsproduct[j] > 0){
              fprintf(OUT,"x%d",rel[i].rhsproduct[j]);
           }
           else if (rel[i].rhsproduct[j] < 0){
              fprintf(OUT,"x%d^(-1)",-rel[i].rhsproduct[j]);
           }
           else{
              fprintf(stderr,"error in put_presentation\n");
              exit(3);
           }

           /* print an intermediate cr for GAP */
           if (j % 20 == 19) fprintf(OUT,"\n");
        }
        if (i != (pres->norelators -1)) fprintf(OUT,",\n");
     }
     fprintf(OUT,"];\n");
     fprintf(OUT,"Size(H);\n");
   }
   else if (strcmp(option,"C") == 0){

      /* the CARAT part */

      /* find the maximal length of a relator */
      for (i=0;i<pres->norelators;i++){
         if (rel[i].lhsnproduct>0){
           fprintf(stderr,"error in put_presentation\n");
           exit(3);
         }
         if (abs(rel[i].rhsnproduct)>max) max = abs(rel[i].rhsnproduct);
      }

      fprintf(OUT,"%dx%d\n",pres->norelators,max);
      for (i=0;i<pres->norelators;i++){
         for (j=0;j<max;j++){
            if (j<rel[i].rhsnproduct)
               fprintf(OUT,"%d ",rel[i].rhsproduct[j]);
            else
               fprintf(OUT,"0 ");
         }
         fprintf(OUT,"\n");
      }
   }

   /* close the file if it was opened */
   if (file != NULL)
      fclose(OUT);

   return;
}  /* put_presentation */

/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ void delete_relator(presentation_TYP *p,int no)
@
@ deletes the relator with number no form the presentation.
@ It will free all records assigned with p->relators[no] if appropriate,
@ and free p->relators[no] as well.
@ The order of the remaining relators won't be changed.
@
@  presentation_TYP *p : The presentation in question.
@  int no              : number of the relator to be deleted
@
@----------------------------------------------------------------------------
@
*****************************************************************************/
void delete_relator(presentation_TYP *p,int no)
{

   int i;

   if (p->norelators <= no){
      fprintf(stderr,"error in delete_relator: no such relator\n");
      exit(3);
   }

   if (p->relators[no].lhsnproduct > 0){
      free(p->relators[no].lhsproduct);
      p->relators[no].lhsproduct = NULL;
      p->relators[no].lhsnproduct = 0;
   }
   if (p->relators[no].rhsnproduct > 0){
      free(p->relators[no].rhsproduct);
      p->relators[no].rhsproduct = NULL;
      p->relators[no].rhsnproduct = 0;
   }

   p->norelators--;

   for (i=no;i<p->norelators;i++){
      memcpy(p->relators+i,p->relators+i+1,sizeof(relator_TYP));
   }

   if (p->norelators == 0){
      free(p->relators);
      p->relators = NULL;
   }

   return;
} /* delete_relator */

/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ void add_generator(presentation_TYP *pres,int no,matrix_TYP *mat,int *word)
@
@ Adds a generator with number 'no' to the presentation pres.
@ If mat!=NULL the function will COPY this matrix to
@ pres->generators->list[no-1]->element.
@ The word is assumed to be an equation satiesfied in the group, ie.
@ new_generator = old_generator[word[1]] * old_generator[word[2]] * ...
@                 old_generator[word[word[0]]]
@
@ presentation_TYP *pres : The presentation in question
@ int no                 : the number of the new generator
@ matrix_TYP *mat        : see above
@ int *word              : a relation satisfied in the group, see above.
@
@ SIDEEFECTS: pres->generators->list will be realloced.
@
@----------------------------------------------------------------------------
@
*****************************************************************************/
void add_generator(presentation_TYP *pres,int no,matrix_TYP *mat,int *word)
{

  int i,
      j;

  pres->generators->firstfree++;

  pres->generators->list = (derived_TYP **) realloc(pres->generators->list,
                          pres->generators->firstfree * sizeof(derived_TYP *)); 
  pres->generators->list[no-1] = (derived_TYP *) malloc(1*sizeof(derived_TYP));

  if (mat == NULL){
     pres->generators->list[no-1]->element = NULL;
  }
  else{
     pres->generators->list[no-1]->element = copy_mat(mat);
  }
  pres->generators->list[no-1]->product = NULL;
  pres->generators->list[no-1]->nproduct = 0;
  pres->generators->list[no-1]->left = 0;
  pres->generators->list[no-1]->right = 0;


  /* fit in the additional relator */
  pres->norelators++;
  if (pres->ext_factor < pres->norelators){
     pres->ext_factor += MIN_SPEICHER;
     pres->relators = (relator_TYP*) realloc(pres->relators,pres->ext_factor
                                        * sizeof(relator_TYP));
  }

  i = pres->norelators-1;
  pres->relators[i].rhsproduct = (int *) malloc((word[0]+1)*sizeof(int));
  for (j=0;j<word[0];j++)
     pres->relators[i].rhsproduct[j] = word[j+1];
  pres->relators[i].rhsproduct[word[0]] = -no;
  pres->relators[i].rhsnproduct = word[0]+1;

  pres->relators[i].lhsnproduct = 0;
  pres->relators[i].lhsproduct = NULL;

  return;
} /* add_generator */

/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ void subs_no_by_word(relator_TYP *relator,int no,int *word,
@                      int **inverse_word)
@
@ substitutes all occurences of 'no' in relator by the word in 'word',
@ and all occurences of '-no' by inverse_word[0].
@ It is possible to set inverse_word[0] = NULL, in this case it will be
@ set to the inverse of 'word' only if it is needed, and the function will
@ return this value in turn.
@
@ relator_TYP *relator: relator in question
@ int no              : the number of the generator to be replaced by word
@                       (or inverse_word if appropriate).
@ int *word           : word of the form x_word[1] * ..... x_word[word[0]].
@ int **inverse_word  : a pointer to a word of the form of word. If
@                       inverse_word[0] == NULL this word will be set to
@                       the inverse of 'word' if in need for this variable.
@
@ NOTE: The function does not assume 'word' and 'inverse_word' to have
@       the same length.
@----------------------------------------------------------------------------
@
*****************************************************************************/
void subs_no_by_word(relator_TYP *relator,int no,int *word,int **inverse_word)
{

   int i,
       j,
       k,
      *tmp,
       number=0,
       inverse_number=0;

  /* count the number of occurences of no ans -no */
  for (i=0;i<relator->rhsnproduct;i++){
    if (relator->rhsproduct[i] == no){
       number++;
    }
    else if (relator->rhsproduct[i] == (-no)){
       inverse_number++;
    }
  }

  /* invert the word if nessesary */
  if ((inverse_number > 0) && inverse_word[0] == NULL){
     inverse_word[0] = (int *) malloc((word[0]+1) * sizeof(int));
     inverse_word[0][0] = word[0];
     for (i=0;i<word[0];i++)
        inverse_word[0][i+1] = -word[word[0]-i];
  }

  /* now we have an upper bound for the length of the resulting relator */
  i = relator->rhsnproduct + word[0]*number + inverse_word[0][0]*inverse_number;
  tmp = (int *) malloc(i * sizeof(int));

  /* subsitute every occurence of 'no' in relator->rhsproduct by the word */
  j=0;
  for (i=0;i<relator->rhsnproduct;i++){
     if (relator->rhsproduct[i] == no){
        for (k=0;k<word[0];k++){
           tmp[j]=word[k+1];
           j++;
        }
     }
     else if(relator->rhsproduct[i] == (-no)){
        for (k=0;k<inverse_word[0][0];k++){
           tmp[j]=inverse_word[0][k+1];
           j++;
        }
     }
     else{
        tmp[j] = relator->rhsproduct[i];
        j++;
     }
  }

  free(relator->rhsproduct);
  relator->rhsnproduct = j;
  relator->rhsproduct = tmp;

  return;
} /* subs_no_by_word */

/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ void delete_generator(presentation_TYP *pres,int no,int *word)
@
@ The function deletes the generator with number 'no' form the presentation.
@ All generators with number > no will have their number decreased by 1.
@ The function will first substitude all occurrences of 'no' (or '-no')
@ in all relators of the presentation by the word 'word' (or its inverse).
@ 'word' is assumed to be a word of the form x_word[1] * ... x_word[word[0]].
@ All relators are afterwards simplified by simplify_relator, and those
@ reducing to the identity word will be cancel. The order of the relators
@ is otherwise left unchanged.
@
@ presentation_TYP *pres: the presentation in question.
@ int no                : the number of the generator to be deleted.
@ int *word             : the word which represent this generator.
@----------------------------------------------------------------------------
@
*****************************************************************************/
void delete_generator(presentation_TYP *pres,int no,int *word)
{
   int i,
       j,
      *inverse_word=NULL;

   /* substitute all occurences of `no' in a relator by the word */
   for (i=0;i<pres->norelators;i++){
      subs_no_by_word(pres->relators+i,no,word,&inverse_word);

      if (simplify_relator(pres->relators+i)){
         delete_relator(pres,i);
         i--; 
      }
   }

   /* all occurences of numbers greater the no have to be decreased by one */
   for (i=0;i<pres->norelators;i++){
      for (j=0;j<pres->relators[i].rhsnproduct;j++){
         if (pres->relators[i].rhsproduct[j]>no){
            pres->relators[i].rhsproduct[j]--;
         }
         else if (pres->relators[i].rhsproduct[j]< (-no)){
            pres->relators[i].rhsproduct[j]++;
         }
      }
   }   

   free(pres->generators->list[no-1]);
   for (i=no;i<pres->generators->firstfree;i++)
      pres->generators->list[i-1] = pres->generators->list[i];
   pres->generators->firstfree--;

   if (inverse_word != NULL) free(inverse_word);

   return;
} /* delete_generator */

/**************************************************************************
@
@--------------------------------------------------------------------------
@
@ static void bubble_sort_relators(presentation_TYP *p)
@
@ Sorts the relators in p according to their length.
@
@ 
@--------------------------------------------------------------------------
@
***************************************************************************/
static void bubble_sort_relators(presentation_TYP *p)
{
   int i;

   relator_TYP *rel=p->relators,
                tmp;

   for (i=1;i<p->norelators;i++){
      if (rel[i-1].rhsnproduct>rel[i].rhsnproduct){
         memcpy(&tmp,rel+i,sizeof(relator_TYP));
         memcpy(rel+i,rel+i-1,sizeof(relator_TYP));
         memcpy(rel+i-1,&tmp,sizeof(relator_TYP));
         i = 0;
      }
   }

   return;
} /*  void bubble_sort_relators(presentation_TYP *p) */

/*************************************************************************
@
@-------------------------------------------------------------------------
@
@ int cancel_subword(relator_TYP *r1,relator_TYP *r2)
@
@ Deletes all subwords r1 in r2. The function assumes that both relators
@ just have a right hand side, and that r1->rhsnproduct <= r2->rhsnproduct.
@ The function will otherwise print an error message to sterr and exit.
@
@-------------------------------------------------------------------------
@
**************************************************************************/
int cancel_subword(relator_TYP *r1,relator_TYP *r2)
{
   int i,
       j;

   if (r1 == r2){
      fprintf(stderr,"error in cancel_subword: arguments have to be\n");
      fprintf(stderr,"different\n");
      exit(3);
   }
   if (r1->rhsnproduct > r2->rhsnproduct){
      fprintf(stderr,"error in cancel_subword: 1st arguments must be\n");
      fprintf(stderr,"smaller than 2nd one.\n");
      exit(3);
   }

   /* look whether r1 is a subword of r2 */
   for (i=0;i<=r2->rhsnproduct-r1->rhsnproduct;i++){
      for (j=0;j<r1->rhsnproduct &&
               (r1->rhsproduct[j] == r2->rhsproduct[i+j]);j++);
      if (j==r1->rhsnproduct){
         r2->rhsnproduct -= r1->rhsnproduct;
         for (j=i;j<r2->rhsnproduct;j++)
           r2->rhsproduct[j] = r2->rhsproduct[j+r1->rhsnproduct];
         return TRUE;
      }
   }

   /* look whether r1^(-1) is a subword of r2 */
   for (i=0;i<=r2->rhsnproduct-r1->rhsnproduct;i++){
      for (j=0;j<r1->rhsnproduct &&
             (-r1->rhsproduct[r1->rhsnproduct-j-1] == r2->rhsproduct[i+j]);j++);
      if (j==r1->rhsnproduct){
         r2->rhsnproduct -= r1->rhsnproduct;
         for (j=i;j<r2->rhsnproduct;j++)
           r2->rhsproduct[j] = r2->rhsproduct[j+r1->rhsnproduct];
         return TRUE;
      }
   }

   return FALSE;

} /* cancel_subword */

/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ void simplify_presentation(presentation_TYP *p)
@
@ Tries to simplify the presentation p. IT WILL NOT CHANGE THE GENERATING
@ SET, and it does only work symbolicly.
@ It repeatedly will call cancel_subword and bubble_sort_relators to
@ archive this goal.
@ NOTE: The function doesn't give a good result in big examples.
@
@----------------------------------------------------------------------------
@
*****************************************************************************/
void simplify_presentation(presentation_TYP *p)
{

   int i,
       j,
       recurse_flag=FALSE;

   relator_TYP *rel=p->relators;

   for (i=0;i<p->norelators;i++){
      if (simplify_relator(rel+i)){
         delete_relator(p,i);
         i--;
      }
   }

   bubble_sort_relators(p);

   for (i=1;i<p->norelators && !recurse_flag;i++){
      for (j=i;j<p->norelators && !recurse_flag;j++){
         recurse_flag = recurse_flag || cancel_subword(rel+i-1,rel+j);
         if (rel[j].rhsnproduct == 0){
            if (rel[j].rhsproduct != NULL) free(rel[j].rhsproduct);
            delete_relator(p,j);
         }
      }
   }

   if (recurse_flag){
      simplify_presentation(p);
   }

   return;
} /*  void simplify_presentation(presentation_TYP *p) */


/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ void tietze(presentation_TYP *pres,int new_anz,int **old_as_new,
@                      int **new_as_old,matrix_TYP **new_mat)
@
@ The function tietze will perform Tietze transformations on the presentation
@ p to change it form one generating set to another.
@ It assumes that all relators are written on the right hand side only.
@ The function will COPY the matrices representing the new generators
@ to the places p->generators->list[??]->element if they are != NULL.
@
@ presentation_TYP *pres: The presentation where the Tietze transformations
@                         will be performed with.
@ int new_anz           : The desired number of generators.
@ int **old_as_new      : for every generator in prer->generators->list
@                         a word w = old_as_new[i] of the form
@                         old_generator_i=new_gen_w[1] *....* new_gen_w[w[0]].
@ int **new_as_old      : for every new generator a word (in the same form
@                         as above), representing it entirely as a product of
@                         generator.
@ matrix_TYP **new_mat  : The matrices representing the new generators.
@                         Note also it is possible to have entries
@                         new_mat[i] == NULL, it is not possible to give
@                         new_mat==NULL. All entries new_mat[i], for
@                         0 <= i < new_anz will be asked for.
@
@----------------------------------------------------------------------------
@
*****************************************************************************/
void tilman_tietze(presentation_TYP *pres,int new_anz,int **old_as_new,
                     int **new_as_old,matrix_TYP **new_mat)
{

   int i,
       j,
       old_anz = pres->generators->firstfree;

   /* add the new generators, easy */
   for (i=0;i<new_anz;i++){
      /* get additional pairs of (generator,relator) */
      add_generator(pres,old_anz+1+i,new_mat[i],new_as_old[i]);
   }

   if (INFO_LEVEL & 4){
      put_presentation(pres,NULL,"G");
   }
/*
   simplify_presentation(pres);
*/

   if (INFO_LEVEL & 4){
      put_presentation(pres,NULL,"G");
   }

   /* remove the new old generators, quite hard */
   for (i=0;i<old_anz;i++){

      /* change the words */
      for (j=1;j<=old_as_new[i][0];j++){
         if (old_as_new[i][j] > 0){
            old_as_new[i][j] += (old_anz - i);
         }
         else{
            old_as_new[i][j] -= (old_anz - i);
         }
      }

      delete_generator(pres,1,old_as_new[i]);

      if (INFO_LEVEL & 4){
         put_presentation(pres,NULL,"G");
      }

      /* rechange the words */
      for (j=1;j<=old_as_new[i][0];j++){
         if (old_as_new[i][j] > 0){
            old_as_new[i][j] -= (old_anz - i);
         }
         else{
            old_as_new[i][j] += (old_anz - i);
         }
      }

   }

   if (INFO_LEVEL & 4){
      put_presentation(pres,NULL,"C");
   }
/*
   simplify_presentation(pres);
*/

   return;
} /* tietze */

/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ presentation_TYP *generate_pres(int no_of_generators)
@
@ allocates all memory to give a presentation on 'no_of_generators'
@  generators.
@
@----------------------------------------------------------------------------
@
*****************************************************************************/
presentation_TYP *generate_pres(int no_of_generators)
{

   presentation_TYP *RES;

   int i;

   RES = (presentation_TYP *) malloc(1 * sizeof(presentation_TYP));

   RES->generators = (derivedsg_TYP *) malloc(1 * sizeof(derivedsg_TYP));
   RES->generators->firstfree = no_of_generators;
   RES->generators->list = (derived_TYP **) malloc(no_of_generators *
                                                   sizeof(derived_TYP *));

   for (i=0;i<no_of_generators;i++){
      RES->generators->list[i] = (derived_TYP *) malloc(1*sizeof(derived_TYP));
      RES->generators->list[i]->element = NULL;
      RES->generators->list[i]->product = NULL;
      RES->generators->list[i]->nproduct = 0;
      RES->generators->list[i]->left = 0;
      RES->generators->list[i]->right = 0;
   }

   RES->relators = 0;
   RES->norelators = 0;
   RES->ext_factor = 0;

   return RES;
} /* generate_pres */
