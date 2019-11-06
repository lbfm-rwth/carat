/***************************************************************************
@
@ --------------------------------------------------------------------------
@
@ FILE: base.h
@
@ some algorithms connected with the schreier-sims algorithm are collected.
@
@ --------------------------------------------------------------------------
@
****************************************************************************/

#include <typedef.h>
#include <longtools.h>
#include <getput.h>
#include <matrix.h>
#include <tools.h>
#include <sort.h>


void free_tree(struct tree *p)
{
   if (p!=NULL){
      free_tree(p->left);
      free_tree(p->right);
      free(p);
   }
}


int hash_mat(struct tree *p,matrix_TYP **v,matrix_TYP *x,int pos)
{
   int cmp;

   if ((p==NULL) || (p->no <0)){
      fprintf(stderr,"fehler in hash_mat\n");
   }

   cmp = mat_comp(v[p->no],x);
   if (cmp < 0){
      if (p->left == NULL){
         if (pos != -1){
            p->left = (struct tree *) calloc(1,sizeof(struct tree));
            p->left->no = pos;
         }
         return -1;
      }
      else{
         return hash_mat(p->left,v,x,pos);
      }
   }
   else if (cmp > 0){
      if (p->right == NULL){
         if (pos != -1){
            p->right = (struct tree *) calloc(1,sizeof(struct tree));
            p->right->no = pos;
         }
         return -1;
      }
      else{
         return hash_mat(p->right,v,x,pos);
      }
   }
   else{
      return p->no;
   }

}

/**************************************************************************
@
@ --------------------------------------------------------------------------
@
@ void init_bahn(bahn *a)
@ initialises the struct bahn, i.e. get's memory for the pointers, and
@ sets all other values to zero.
@
@ --------------------------------------------------------------------------
@
***************************************************************************/
extern void init_bahn(bahn *a)
{
   a->length = 0;
   a->speicher = MIN_SPEICHER;
   a->orbit = (matrix_TYP **) calloc(MIN_SPEICHER , sizeof(matrix_TYP *));
   a->representatives = (matrix_TYP **) calloc(MIN_SPEICHER
                                         , sizeof(matrix_TYP *));
   a->rep_invs = (matrix_TYP **) calloc(MIN_SPEICHER
                                         , sizeof(matrix_TYP *));
   a->words = (int **) calloc(MIN_SPEICHER, sizeof(int *));
   a->generators = (matrix_TYP **) calloc(MIN_SPEICHER
                                         , sizeof(matrix_TYP *));
   a->hash = (struct tree *) malloc(1*sizeof(struct tree));
   a->hash->left = NULL;
   a->hash->right = NULL;
   a->hash->no = -1;

}


/**************************************************************************
@
@ --------------------------------------------------------------------------
@
@ extern void extend_bahn(bahn **a)
@
@  bahn **a
@
@  reallocated the pointers a[0]->orbit, a[0]->generators and
@  a[0]->representatives by the amount MIN_SPEICHER, and increases
@  a[0]->speicher by this ammount.
@
@ --------------------------------------------------------------------------
@
***************************************************************************/
extern void extend_bahn(bahn **a)
{
   if (a[0]->speicher == 0){
      printf("error in extend_bahn\n");
      exit(3);
   }
   a[0]->speicher = a[0]->speicher+MIN_SPEICHER;
   a[0]->orbit = (matrix_TYP **) realloc(a[0]->orbit,
                 sizeof(matrix_TYP *) * a[0]->speicher);
   a[0]->representatives = (matrix_TYP **) realloc(a[0]->representatives,
                 sizeof(matrix_TYP *) * a[0]->speicher);
   a[0]->rep_invs = (matrix_TYP **) realloc(a[0]->rep_invs,
                 sizeof(matrix_TYP *) * a[0]->speicher);
   a[0]->words = (int **) realloc(a[0]->words,sizeof(int*));
   a[0]->generators = (matrix_TYP **) realloc(a[0]->generators,
                 sizeof(matrix_TYP *) * a[0]->speicher);
}

/***************************************************************************
@
@ --------------------------------------------------------------------------
@
@ void free_bahn(bahn *a)
@
@ --------------------------------------------------------------------------
@
***************************************************************************/
extern void free_bahn(bahn *a)
{
   int i;

   for (i=0;i<a->length;i++){
     if (a->orbit[i] != NULL){
        free_mat(a->orbit[i]);
     }
     if (a->representatives[i] != NULL){
        free_mat(a->representatives[i]);
     }
     if (a->rep_invs[i] != NULL){
        free_mat(a->rep_invs[i]);
     }
     if (a->words[i] != NULL){
        free(a->words[i]);
     }
   }
   a->length =0;
   a->gen_no = 0;
   a->speicher = 0;
   free(a->orbit);
   free(a->representatives);
   free(a->rep_invs);
   free(a->words);
   free(a->generators);
   a->orbit = NULL;
   a->representatives = NULL;
   a->generators = NULL;

   free_tree(a->hash);
   a->hash = NULL;
}

static int position(matrix_TYP **a,matrix_TYP *x,int n)
/* returns the first index i<n such that a[i] == x, and
   -1 if none exists */
{
  int i=0;

  while (i<n){
    if (mat_comp(a[i],x) == 0){
       return i;
    }
    i++;
  }
  return -1;
}

static void inv_word(int *w,
                     int *w_to_inv)
{
   int i,
       j;

   for (i=1,j=w_to_inv[0];i<=w_to_inv[0];i++,j--)
      w[i] = -w_to_inv[j];

   w[0] = w_to_inv[0];

   return;
} /* inv_word */

static int *get_word_to_generator(bahn *A,
                                  matrix_TYP *h)
{

  int i,
     *w;

  for (i=0;i<A->length && mat_comp(A->representatives[i],h) ; i++);

  w = (int *) malloc((A->words[i][0]+1)*sizeof(int));
  memcpy(w,A->words[i],(A->words[i][0]+1)*sizeof(int));

  return w;

} /* get_word_to_generator */



/***************************************************************************
@
@ --------------------------------------------------------------------------
@
@  int einordnen(bahn** erg,
@                matrix_TYP *h,
@                matrix_TYP *new_vec,
@                int *w,
@                int l,
@                int flag)
@
@  bahn **erg          :  the set of base/strong generators for a finite
@                         group G. it's the ONLY variable to be changed.
@  matrix_TYP *h       :
@  matrix_TYP *new_vec :  this value is not asked for in the case
@                         described here, so any variable with the right
@                         type will do.
@  int l               :  -1          ALLWAYS
@  int flag            :  TRUE        ALLWAYS!!!!!
@
@ Inserts a new generator h into the list of strong generators for
@ any FINITE group G.
@ Therefore erg has to be the result of a call of strong_generators(..)
@ for some bravais_TYP *G.
@ The function calls itself recursively, and therefore the convention to
@ use it is somewhat cryptical. In all cases for the 'end user' l has to
@ be -1, and flag has to be TRUE (no garanties otherwise).
@ THE FUNCTION DOES NOT COPY THE VARIABLE h INTO erg[0]->generators,
@ (IT JUST STORES THE POINTER THERE), THEREFORE THE USER IS ASKED TO
@ ASSURE THAT h STAY'S IN THE SAME PLACE WHILE USING erg !!!!!!!!!!
@ The function does not check whether the new generator h is really needed,
@ so better check is out before calling this function by a call of is_element!
@ --------------------------------------------------------------------------
@
***************************************************************************/
static int einordnen(bahn** erg,
                     matrix_TYP *h,
                     matrix_TYP *new_vec,
                     int *w,
                     int l,
                     int flag)
{
   int tmp = erg[0]->representatives[0]->cols,
       tmp2,
       i,
      *wn;

   matrix_TYP  *hinv,
               *nv,
               *nh;

   if (INFO_LEVEL & 2){
      fprintf(stderr,"entering einordnen with l = %d\n",l);
   }

   if (l>=tmp){
      return FALSE;
   }
   else{

      if (l!=(-1)){
         tmp2 = hash_mat(erg[l]->hash,erg[l]->orbit,new_vec,erg[l]->length);
      }
      else{
         tmp2 = (-1);
      }

      if (tmp2!= (-1)){
         /* we got an element of the stabilizer, */
         /* calculate it then!                   */
         free_mat(new_vec);

         /* in the last stage we won't need an element of
            the stabilizer anymore, so we won't calculate it */
         if (((l+1)<tmp) && (mat_comp(h,erg[l]->representatives[tmp2])!=0)){
            /* hinv = mat_inv(h);
            free_mat(h);
            mat_muleq(hinv,erg[l]->representatives[tmp2]);
            h = hinv;
            new_vec = mat_mul(h,erg[l+1]->orbit[0]); */
            hinv = mat_mul(erg[l]->rep_invs[tmp2],h);
            free_mat(h);
            h = hinv;
            new_vec = mat_mul(h,erg[l+1]->orbit[0]);

            if (w){
	      wn = (int *) malloc((erg[l]->words[tmp2][0]+w[0]+1)
                             *sizeof(int));
              inv_word(wn,erg[l]->words[tmp2]);
              memcpy(wn+wn[0]+1,w+1,w[0]*sizeof(int));
              wn[0] = wn[0] + w[0];
              free(w);
            }
            else{
               wn = NULL;
            }
         }
         else{
            /* only clean up memory */
            free_mat(h);
            if (w) free(w);
            return FALSE;
         }

         /* enter the next stage */
         flag = einordnen(erg,h,new_vec,wn,++l,TRUE);
      }
      else{
         if (l != (-1)){
            if (erg[l]->speicher<=erg[l]->length){
               extend_bahn(erg+l);

               /* this is so seldom used that we might as well check wheter
                  the group is infinite, find a solution later !! */
            }
            erg[l]->orbit[erg[l]->length]=new_vec;
            erg[l]->representatives[erg[l]->length]=h;
            /* inserted to reduce the number of mat_inv s */
            erg[l]->rep_invs[erg[l]->length] = mat_inv(h);

            if (w){
               /* deal with the words */
	       erg[l]->words[erg[l]->length] 
                     = (int *) malloc((w[0]+1) * sizeof(int));
               memcpy(erg[l]->words[erg[l]->length],w,(w[0]+1) * sizeof(int));
               free(w);
	    }
            erg[l]->length++;

            if (INFO_LEVEL & 1){
               fprintf(stderr,"l = %d, lenght in this stage %d\n"
                           ,l,erg[l]->length);
            }
         }

         if (flag){

            if (l == (-1)) l =0;

            /* firstly put this element on the generator list */
            erg[l]->generators[erg[l]->gen_no] = h;
            erg[l]->gen_no++;

         }

      }
   }

   if (flag){
      /* do a new orbit-calculation */
      /* apply all generators to the previous orbit elements */
      for(tmp=0;tmp<erg[l]->length;tmp++){
         for (tmp2=l;tmp2<erg[0]->representatives[0]->cols;tmp2++){
            for (i=erg[tmp2]->gen_no-1;i>=0;i--){
               h = erg[tmp2]->generators[i];
               nv = mat_mul(h,erg[l]->orbit[tmp]);
               nh = mat_mul(h,erg[l]->representatives[tmp]);
               if (w){
                  wn = get_word_to_generator(erg[tmp2],h);
                  wn = (int *) realloc(wn,
                             (wn[0]+erg[l]->words[tmp][0]+1) * sizeof(int));
                  memcpy(wn+wn[0]+1,erg[l]->words[tmp]+1,
                             erg[l]->words[tmp][0] * sizeof(int));
                  wn[0] = wn[0] + erg[l]->words[tmp][0];
               }
               else{
                  wn = NULL;
               }
               einordnen(erg,nh,nv,wn,l,FALSE);
            }
         }
      }
   }

   return flag;
} /* einordnen(..........) */

/***********************************************************************
@
@-----------------------------------------------------------------------
@
@ bahn **strong_generators(matrix_TYP **base,bravais_TYP *U)
@
@  matrix_TYP **base:
@  bravais_TYP *U   :
@
@  Calculates a strong generating set according to the sims algorithm
@  for the FINITE group U and the given base.
@  This function initialises the result according to the conventions
@  made in this file. It doesn't change neither **base nor *U.
@  for U only the entries U->dim, U->gen and U->gen_no are used.
@-----------------------------------------------------------------------
@
************************************************************************/
bahn **strong_generators(matrix_TYP **base,bravais_TYP *U,int OPT)
{
   bahn **erg;

   int i,
       j,
       k,
       m,
       dim = U->dim,
      *w;

   matrix_TYP **gen=U->gen,
               *new_vec,
               *g,
               *h;


   /* geting the memory for the result */
   erg = (bahn **) calloc(dim,sizeof(bahn *));
   for (i=0;i<dim;i++){
      erg[i] = (bahn *) calloc(1,sizeof(bahn));
      init_bahn(erg[i]);
   }

   for (i=0;i<dim;i++){
      erg[i]->representatives[0] = init_mat(U->dim,U->dim,"1");
      erg[i]->rep_invs[0] = init_mat(U->dim,U->dim,"1");
      erg[i]->words[0] = (int *) malloc(sizeof(int));
      erg[i]->words[0][0] = 0;
      erg[i]->orbit[0] = copy_mat(base[i]);
      erg[i]->length = 1;
      erg[i]->hash->no = 0;
   }

   /* doing the orbit-calculation */
   for (j=0;j<erg[0]->length;j++){
     for (k=0;k<U->gen_no;k++){

        if (INFO_LEVEL & 4){
        /* printing all results for debugging */
        for (i=0;i<dim;i++){
           for (m=0;m<erg[i]->length;m++){
              printf("i=%d m=%d\n",i,m);
              if (erg[i]->representatives[m]->rows != dim ||
                  erg[i]->orbit[m]->rows != dim){
                printf("fehler in strong_generators\n");
                exit(3);
              }
              put_mat(erg[i]->representatives[m],NULL,"represe",2);
              put_mat(erg[i]->orbit[m],NULL,"orbit",2);
           }
        }
        }

        /* decide which group element will act */
        g=gen[k];

        if (INFO_LEVEL & 4 ){
           put_mat(g,NULL,"angewandter erzeuger",2);
        }

        h = mat_mul(g,erg[0]->representatives[j]);
        new_vec = mat_mul(g,erg[0]->orbit[j]);


        if (OPT){
           w = (int *) malloc((erg[0]->words[j][0] + 2) * sizeof(int));
           w[0] = erg[0]->words[j][0] + 1;
           w[1] = k+1;
           memcpy(w+2,erg[0]->words[j]+1,erg[0]->words[j][0] * sizeof(int));
        }
        else{
           w = NULL;
	}
        einordnen(erg,h,new_vec,w,0,FALSE);

      }
    }

   return erg;
} /* strong_generators */

extern matrix_TYP **get_base(bravais_TYP *U)
/* will only return a standart basis */
{
  matrix_TYP **erg;

  int dim=U->dim,
      i;

  erg = (matrix_TYP **) malloc(dim * sizeof(matrix_TYP *));

  for (i=0;i<dim;i++){
    erg[i] = init_mat(dim,1,"i");
    erg[i]->array.SZ[i][0] = 1;
  }

  return erg;
} /* get_base */

/*************************************************************************
@
@ --------------------------------------------------------------------------
@
@ int is_element(matrix_TYP *x,bravais_TYP *G,bahn **strong,int **w)
@
@      matrix_TYP *x:  the matrix in question
@      bravais_TYP *G: actually the only thing asked for is G->dim!
@      bahn **strong:  the result of a call of strong_generators(..)
@                      for the group G.
@
@  Decides whether the matrix represented by *x is contained in the
@  group described by **strong and *G respectively.
@  strong has to be the result of a call of strong_generators(..)
@  for the group G.
@  The function returns TRUE if x in G, and FALSE otherwise.
@  None of the given variables is changed.
@
@ --------------------------------------------------------------------------
@
***************************************************************************/
extern int is_element(matrix_TYP *x,bravais_TYP *G,bahn **strong,int **w)
{
  int erg=TRUE,
      dim = G->dim,
      i=0,
      pos,
     *w2;

  matrix_TYP *g_neu,
             *bild,
             *h;

  g_neu = copy_mat(x);

  while (erg && (i<dim)){

     /* map the basepoint by the element which got to be in
        the pointwise stabilizer of the basepoints b_0,..b_{i-1} */
     bild = mat_mul(g_neu,strong[i]->orbit[0]);

     /* test whether this image is an image also produced by G */
     pos = hash_mat(strong[i]->hash,strong[i]->orbit,bild,-1);

     if (pos == (-1)){
       erg = FALSE;
       free_mat(bild);
     }
     else{
       /* h = mat_inv(g_neu);
       free_mat(g_neu);
       free_mat(bild);
       g_neu = mat_mul(h,strong[i]->representatives[pos]);
       free_mat(h); */
       free_mat(bild);
       h = mat_mul(strong[i]->rep_invs[pos],g_neu);
       free_mat(g_neu);
       g_neu = h;

       if (w){
	 if (w[0]){
	   w2 = (int *) calloc(w[0][0]+strong[i]->words[pos][0]+1, 
			       sizeof(int));
           w2[0] = w[0][0]+strong[i]->words[pos][0];
           memcpy(w2+1,w[0]+1,w[0][0] * sizeof(int));
           memcpy(w2+w[0][0]+1,(strong[i]->words[pos])+1,
		  strong[i]->words[pos][0]* sizeof(int));
           free(w[0]);
	 }
         else{
           w2 = (int *) calloc(strong[i]->words[pos][0]+1,sizeof(int));
           w2[0] = strong[i]->words[pos][0];
           memcpy(w2+1,strong[i]->words[pos]+1,strong[i]->words[pos][0]*
		  sizeof(int));
	 }
         w[0] = w2;
       }
       i++;

     }
  }

  free_mat(g_neu);

  return erg;
} /* is_element */

/**************************************************************************
@
@ --------------------------------------------------------------------------
@
@  int size(bahn **a)
@
@  bahn **a:   the result of a call strong_generators(...) for some
@              bravais_TYP *G.
@
@ Under this condition, the function will return the order of the group G.
@ Actually calculated is the product of a[i]->length in some ranges.
@
@ --------------------------------------------------------------------------
@
**************************************************************************/
extern int size(bahn **a)
{
  int i,
      size = 1;

  for (i=0;i<a[0]->representatives[0]->cols;i++){
     size = size*a[i]->length;
  }

  return size;
}

/**************************************************************************
@
@ --------------------------------------------------------------------------
@
@  matrix_TYP **normalizer_in_N(bravais_TYP *U,bravais_TYP *N,int *anz,
@                                    int finite)
@
@  calculates a generating set for the normalizer of the FINITE group
@  U in the group N, and N is interpreted as generated by N->gen, N->cen and
@  N->normal. The assumption is, that U^N is FINITE.
@     bravais_TYP *U
@     bravais_TYP *N
@     int *anz
@     int finite
@  It returns the elements of the normalizer via return, and the
@  number of these via anz[0].
@  no global variables nor U or N are changed.
@  If finite == true the program assumes that the calculated normalizer
@  is finite, and tries to get hold of a set of base/strong generators
@  for it. THIS WILL WORK IFF THE CALCULATED NORMALIZER IS INDEED FINITE.
@  The program will get the record U->order and U->divisors!
@
@ --------------------------------------------------------------------------
@
***************************************************************************/
extern matrix_TYP **normalizer_in_N(bravais_TYP *U,bravais_TYP *N,int *anz,
                                    int finite)
{
   matrix_TYP **erz_N,
              **erg,
              **orbit,
               *test_subgroup,
              **base,
               *tmp,
               *tmp2,
               *tmp3,
               *stab_element;

   bravais_TYP BN;

   bahn **strong,
        **strong_N = 0;

   int i,
       j,
       k,
       l,
       flag,
       gefunden,
       gen_no,
       dim = U->dim,
       speicher,
       erg_speicher;

   /* getting the standart basis for the vectorspace */
   base = get_base(U);
   strong = strong_generators(base,U,FALSE);

   /* set the order of U */
   U->order = size(strong);
   factorize_new(U->order,U->divisors);
 
   if (INFO_LEVEL & 4){
      printf("calculated the base for U\n");
   }

   /* we are going to need the inverse of each generator anyway
      for conjugation, so calculate them now and stick them
      in the right position */
   gen_no = N->gen_no + N->normal_no + N->cen_no;

   erz_N = (matrix_TYP **) malloc(gen_no * sizeof(matrix_TYP *));

   for (i=0;i<gen_no;i++){
      if (i<N->gen_no){
         erz_N[i] = copy_mat(N->gen[i]);
      }
      else if(i<N->gen_no + N->cen_no){
         erz_N[i] = copy_mat(N->cen[i - N->gen_no]);
      }
      else {
         erz_N[i] = copy_mat(N->normal[i-N->gen_no-N->cen_no]);
      }
   }

   /* the subgroups will be represented by the element that
      conjugates them  */
   orbit = (matrix_TYP **) malloc(MIN_SPEICHER * sizeof(matrix_TYP *));
   speicher = MIN_SPEICHER;
   orbit[0] = init_mat(dim,dim,"1");
   gefunden = 1;

   /* getting memory for the result */
   erg = (matrix_TYP **) malloc(MIN_SPEICHER * sizeof(matrix_TYP *));
   erg_speicher = MIN_SPEICHER;
   anz[0] = 0;

   /* strong orbit calculation yields the normalizer */
   for (i=0;i<gefunden;i++){

      if (INFO_LEVEL & 1){
         fprintf(stderr,"found %d conjugated subgroups,\n",gefunden);
         fprintf(stderr,"dealt with %d of these.\n",i);
      }

      for (j=0;j<(gen_no);j++){

         /* we only have do deal with the generators, because, we
            assume that the orbits are finite                     */
         test_subgroup = mat_mul(erz_N[j],orbit[i]);

         /* now look, whether we already got this subgroup */
         /* it's a fiddly part, because I can't think of a */
         /* ordering for these subgroups                   */
         k=0;
         flag = FALSE;
         tmp = mat_inv(test_subgroup);

         while (!flag && (k<gefunden)){

            if (INFO_LEVEL & 4){
              put_mat(tmp,NULL,"tmp",2);
              put_mat(orbit[k],NULL,"orbit[k]",2);
              printf("k= %d\n",k);
            }
            stab_element = mat_mul(tmp,orbit[k]);

            if (INFO_LEVEL & 4){
               put_mat(stab_element,NULL,NULL,2);
            }

            /* calculate the conjugates of the generators of U
               and decide whether they are in U */
            l=0;
            while (l<U->gen_no){
               if (INFO_LEVEL & 4){
                  put_mat(stab_element,NULL,"stab_el",2);
               }
               tmp2 = mat_inv(stab_element);
               tmp3 = mat_mul(tmp2,U->gen[l]);
               free_mat(tmp2);
               tmp2 = mat_mul(tmp3,stab_element);
               free_mat(tmp3);

               if (INFO_LEVEL & 4){
                  printf("l= %d\n",l);
               }

               if (!is_element(tmp2,U,strong,NULL)){
                 l = U->gen_no + 1;
               }
               l++;
               free_mat(tmp2);
            }

            if (l==U->gen_no){
               flag=TRUE;
            }
            else{
               free_mat(stab_element);
               stab_element = NULL;
            }
            k++;

         }/* while (!flag .. */

         free_mat(tmp);

         /* so now we see whether the element was really in the
            stabilizer */
         if (flag){
            free_mat(test_subgroup);

            /* we might already got this element of the normalizer,
               so better check it out */
            if ((position(erg,stab_element,anz[0])> (-1)) &&
                 (stab_element != NULL)){
               free_mat(stab_element);
            }
            else if (finite){
               /* if we assume the normalizer to be finite, we will
                  have an even better test for a new element */
               if (anz[0] == 0){
                  erg[0] = stab_element;
                  anz[0]++;
                  BN.gen = erg;
                  BN.gen_no = 1;
                  BN.dim = dim;
                  strong_N = strong_generators(base,&BN,FALSE);
                  strong_N[0]->generators[0] = stab_element;
                  strong_N[0]->gen_no = 1;
               }
               else{
                  /* we might already got this element in our stabilizer */
                  if (is_element(stab_element,&BN,strong_N,NULL)){
                     free_mat(stab_element);
                     stab_element = NULL;
                  }
                  else{
                     /* this element is new */
                     anz[0]++;
                     if (erg_speicher <= anz[0]){
                        erg_speicher = erg_speicher+MIN_SPEICHER;
                        erg = (matrix_TYP **) realloc(erg,
                                erg_speicher*sizeof(matrix_TYP *));
                     }
                     erg[anz[0]-1] = stab_element;
                     einordnen(strong_N,stab_element,
                               stab_element,NULL,-1,TRUE);
                  }
               }
            }
            else{
               anz[0]++;
               if (erg_speicher<=anz[0]){
                  erg_speicher = erg_speicher+MIN_SPEICHER;
                  erg = (matrix_TYP **) realloc(erg,
                                erg_speicher*sizeof(matrix_TYP *));
               }
               erg[anz[0]-1] = stab_element;
            }
         }
         else{
            if (INFO_LEVEL & 4){
               put_mat(stab_element,NULL,NULL,2);
            }

            if (stab_element != NULL) free_mat(stab_element);

            gefunden++;
 
            if (speicher<=gefunden){
               speicher = speicher+MIN_SPEICHER;
               orbit = (matrix_TYP **) realloc(orbit,
                                 speicher*sizeof(matrix_TYP *));
            }
            orbit[gefunden-1] = test_subgroup;
         }

      }/* for (j=0 .. */
   }/* for (i=0 .., ie the orbit calculation */

   /* cleaning up the memory */
   for (i=0;i<dim;i++){
      free_mat(base[i]);
      free_bahn(strong[i]);
      free(strong[i]);
      if (finite){
         free_bahn(strong_N[i]);
         free(strong_N[i]);
      }
   }
   if ((INFO_LEVEL & 1) && finite){
      fprintf(stderr,"Order of the normalizer: %d\n",size(strong_N));
   }
   free(base);
   free(strong);
   if (finite) free(strong_N);

   for (i=0;i<gen_no;i++){
      free_mat(erz_N[i]);
   }
   free(erz_N);

   for (i=0;i<gefunden;i++){
      free_mat(orbit[i]);
   }
   free(orbit);

   return erg;
} /* normalizer_in_N */
 
/**************************************************************************
@
@--------------------------------------------------------------------------
@
@ int red_gen(bravais_TYP *G,matrix_TYP **base,bahn ***strong,int i)
@
@ Tries to reduce the number of generators for the group G, and returns
@ |G|.
@
@   bravais_TYP *G:     the group in question
@   matrix_TYP **base:  a basis for the vectorspace on which G is acting
@   bahn ***strong:     after the call of this function strong[0] will
@                       contain a set of base/strong generating set for
@                       the group G.
@                       There are two cases for the call of this function
@                       1st: strong[0] == NULL
@                            the pfunction reduces the number of generators
@                            for G.
@                       2nd: strong[0] != NULL
@                            the function assummes strong[0] to contain a set
@                            of base/strong generators for the group generated
@                            by the first i generators of G, and kills
@                            those of the other generators which does not
@                            give a bigger group.
@   int i:              see above
@
@   SIDEEFFECTS:        the entries of G->gen may be changed, and they
@                       may be shuffled. G->gen_no may be changed to
@                       the number of relevant matrices.
@                       strong[0] returns a set of base/strong generators
@                       for G, NO OTHER GLOBAL VARIABLE IS CHANGED.
@
@--------------------------------------------------------------------------
@
***************************************************************************/
extern int red_gen(bravais_TYP *G,matrix_TYP **base,bahn ***strong,int i)
{
   int G_gen_no= G->gen_no;

   if (strong[0] == NULL){
      G->gen_no = 1;
      i=1;
      strong[0] = strong_generators(base,G,FALSE);
      strong[0][0]->gen_no = 1;
      strong[0][0]->generators[0] = G->gen[0];
   }

   /* this will free all irrelevant generators */
   while (i<G_gen_no){
      if (!is_element(G->gen[i],G,strong[0],NULL)){
         einordnen(strong[0],G->gen[i],G->gen[i],NULL,-1,TRUE);
      }
      else{
         free_mat(G->gen[i]);
         G->gen[i] = NULL;
      }
      i++;
   }

   G->gen_no = G_gen_no;
   for (i=0;i<G->gen_no;i++){
      if (G->gen[i] == NULL){
         G->gen[i] = G->gen[G->gen_no-1];
         G->gen_no--;
         i--;
      }
   }

   /* just for sanity */
   G->gen = (matrix_TYP **) realloc(G->gen, G->gen_no * sizeof(matrix_TYP *));

   return size(strong[0]);
}
