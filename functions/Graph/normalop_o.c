/****************************************************************************
@
@ ----------------------------------------------------------------------------
@
@ FILE: normalop.c
@
@ ----------------------------------------------------------------------------
@
******************************************************************************/

#include <typedef.h>
#include <getput.h>
#include <matrix.h>
#include <base.h>
#include <gmp.h>
#include <zass.h>
#include <longtools.h>
#include <orbit.h>
#include <sort.h>
#include <bravais.h>
#include <contrib.h>
#include <graph.h>
#include <datei.h>

#define TWOTO21 2097152

extern int INFO_LEVEL;

static void mod(int *a,int b)
{
   a[0] = a[0] % b;
   if (a[0]<0) a[0] += b;
}

static void valuation_here(matrix_TYP *x,matrix_TYP *D,MP_INT *val)
{
   int first,
       last,
       i;

   MP_INT prod,
          tmp;

   mpz_init_set_si(&prod,1);
   mpz_set_si(val,0);
   mpz_init(&tmp);

   /* set first and last */
   for (first = 0;first<D->cols && D->array.SZ[first][first] == 1;first++);
   for (last = 0;last<D->cols && D->array.SZ[last][last] != 0;last++);

   for (i=first;i<last;i++){
     mod(x->array.SZ[i-first],D->array.SZ[i][i]);
     mpz_mul_ui(&tmp,&prod,(unsigned long ) x->array.SZ[i-first][0]);
     mpz_add(val,val,&tmp);
     mpz_mul_ui(&prod,&prod,(unsigned long) D->array.SZ[i][i]);
   }

   mpz_clear(&prod);
   mpz_clear(&tmp);

   return;
} /* valuation(.....) */

static int hash(struct tree *p,MP_INT *valuations,MP_INT *new_val,int pos)
{

   int cmp;

   if ((p==NULL) || (p->no <0)){
      fprintf(stderr,"fehler in hash\n");
   }

   cmp = mpz_cmp(new_val,valuations+p->no);

   if (cmp < 0){
      if (p->left == NULL){
         if (pos != -1){
            p->left = (struct tree *) calloc(1,sizeof(struct tree));
            p->left->no = pos;
         }
         return -1;
      }
      else{
         return hash(p->left,valuations,new_val,pos);
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
         return hash(p->right,valuations,new_val,pos);
      }
   }
   else{
      return p->no;
   }

} /* static int hash(...) */

static int smallest(struct tree *p)
{
   if (p==NULL){
      printf("error in smallest\n");
      exit(3);
   }

   if (p->left == NULL){
     return p->no;
   }
   else{
     return smallest(p->left);
   }

} /* static int smallest(...) */

/**************************************************************************
@
@ -------------------------------------------------------------------------
@
@ matrix_TYP *orbit_rep(matrix_TYP *x,
@                       matrix_TYP **N,
@                       int nanz,
@                       matrix_TYP *D,
@                       int option,
@                       char *B,
@                       MP_INT *l,
@                       int *anz,
@                       int **word,
@                       int word_flag,
@                       int ***WORDS,
@                       int *NUMBER_OF_WORDS)
@
@
@ option:       an integer controling the behaviours of the function.
@               option = 0: return the length of the orbit via anz[0] and
@                           the smallest representative via return.
@               option = 1: return straight if we found a smaller
@                           representative, and return this one.
@                           anz[0] has no defined meaning in this case.
@ char *B:
@
@ l:            the valuation of the representative returned is via
@               this pointer.
@
@ int **word:   needed only for word_flag==TRUE:
@               returns a list of integers word with the following property:
@               N[word[0][1]] * N[word[0][2]] * ... *N[word[0][word[0][0]]
@               is a word representing the smallest representative in the
@               orbit.
@ int word_flag: drives the behaviour of the function:
@                for word_flag == TRUE it calculates such a word described
@                above, otherwise it will not change word[0] (hopefully).
@ int ***WORDS:
@ int *NUMBER_OF_WORDS:
@
@ Neither the matrices given nor any global variable is changed.
@
@ -------------------------------------------------------------------------
@
***************************************************************************/
static matrix_TYP *static_orbit_rep(matrix_TYP *x,
                             matrix_TYP **N,
                             int nanz,
                             matrix_TYP *D,
                             int option,
                             char *B,
                             MP_INT *l,
                             int *anz,
                             int **word,
                             int word_flag,
                             int ***WORDS,
                             int *NUMBER_OF_WORDS)
{
   int first,
       last,
       i,
       j,
       k,
       h,
     **orb_words,
       speicher;      /* the number of allocated memory for valuations
                          and orbit, and orb_words */
   MP_INT *valuations,
           new_val;        /* the valuation of the newly calculated vector */

   matrix_TYP *erg,
             **orbit,
              *tmp;

   struct tree *p;

   /* set first and last */
   for (first = 0;first<D->cols && D->array.SZ[first][first] == 1;first++);
   for (last = 0;last<D->cols && D->array.SZ[last][last] != 0;last++);

   /* get memory */
   p = (struct tree *) calloc(1,sizeof(struct tree));
   speicher = MIN_SPEICHER;
   valuations = (MP_INT *) malloc(speicher * sizeof(MP_INT));
   for (i=0;i<speicher;i++) mpz_init(&valuations[i]);
   mpz_init(&new_val);
   orbit = (matrix_TYP **) malloc(speicher * sizeof(matrix_TYP*));
   orbit[0] = copy_mat(x);
   if (word_flag){
      orb_words = (int **) malloc(speicher * sizeof(int ));
      orb_words[0] = (int *) malloc(sizeof(int));
      orb_words[0][0] = 0;
   }
   valuation_here(x,D,&valuations[0]);

   /* anz[0] is the length of the orbit so far */
   anz[0] = 1;
   erg = NULL;

   for (i=0;i<anz[0] && erg == NULL;i++){
      for (j=0;j<nanz && erg == NULL;j++){
         tmp = mat_mul(N[j],orbit[i]);

         /* normalize tmp modulo the diagonal entries of D */
         for (k=first;k<last;k++)
           mod(tmp->array.SZ[k-first],D->array.SZ[k][k]);

         valuation_here(tmp,D,&new_val);

         h = hash(p,valuations,&new_val,anz[0]);

         if (h != (-1)){
            /* the element is not new */
            free_mat(tmp);

            /* inserted this case to handle the stabilizer of a cozycle
               as well: tilman 15.03. */
            if (WORDS){

               /* we got a new generator of the stabilizer of the cozycle */
               WORDS[0][NUMBER_OF_WORDS[0]] = (int *) calloc(
                             orb_words[h][0]+orb_words[i][0]+2, sizeof(int));
               WORDS[0][NUMBER_OF_WORDS[0]][0] =
                              orb_words[h][0]+orb_words[i][0]+1;
               for (k=1;k<=orb_words[h][0];k++){
                  WORDS[0][NUMBER_OF_WORDS[0]][k] =
                           -orb_words[h][orb_words[h][0] - k + 1] - 1;
               }
               WORDS[0][NUMBER_OF_WORDS[0]][orb_words[h][0]+1] = j+1;
               for (k=1;k<=orb_words[i][0];k++){
                  WORDS[0][NUMBER_OF_WORDS[0]][orb_words[h][0]+1+k] =
                           orb_words[i][k] + 1;
               }

               NUMBER_OF_WORDS[0]++;
               if (NUMBER_OF_WORDS[0] % MIN_SPEICHER)
                  WORDS[0] = (int **) realloc( WORDS[0] ,
                               (NUMBER_OF_WORDS[0]+MIN_SPEICHER) * sizeof(int));
            }
         }
         else{
            /* the element is new */

            /* get more memory is nessesary */
            if (anz[0] >= speicher){
               speicher += MIN_SPEICHER;
               valuations = (MP_INT *) realloc(valuations,
                                       speicher * sizeof(MP_INT));
               for (k=anz[0];k<speicher;k++) mpz_init(&valuations[k]);
               orbit = (matrix_TYP **) realloc(orbit,
                                    speicher * sizeof(matrix_TYP *));
               if (word_flag){
                  orb_words = realloc(orb_words,speicher * sizeof(int));
               }
            }

            mpz_set(valuations+anz[0],&new_val);
            orbit[anz[0]] = tmp;
            anz[0]++;

            if (word_flag){
               orb_words[anz[0]-1] = (int *) malloc( (orb_words[i][0] + 2)
                                                    * sizeof(int));
               memcpy(orb_words[anz[0]-1]+2,orb_words[i]+1,(orb_words[i][0])
                                                       * sizeof(int));
               orb_words[anz[0]-1][0] = orb_words[i][0] + 1;
               orb_words[anz[0]-1][1] = j;
            }

            /* here is the place to check whether the option allows us
               to exit straigth if tmp is small enough */
            if ((option == 1) && (mpz_cmp(&new_val,&valuations[0])<0)){
               erg = tmp;
               mpz_set(l,&new_val);
               anz[0]--;
            }

            /* mark this element got in the element list B */
            if (B!= NULL) B[mpz_get_ui(&new_val)] = TRUE;
         }
      }
   }

   /* now get the smallest representative if nessesary */
   if (option == 0 || (option == 1 && erg == NULL)){
      k = smallest(p);
      erg = copy_mat(orbit[k]);
      mpz_set(l,&valuations[k]);
      if (word_flag){
         word[0] = orb_words[k];
         /* just malloced to be able to free it */
         orb_words[k] = (int *) malloc(1*sizeof(int));
      }
   }

   /* cleaning up */
   for (i=0;i<anz[0];i++) free_mat(orbit[i]);
   free(orbit);
   for (i=0;i<speicher;i++) mpz_clear(valuations+i);
   free(valuations);
   free_tree(p);
   mpz_clear(&new_val);
   if (word_flag){
      for (i=0;i<anz[0];i++) free(orb_words[i]);
      free(orb_words);
   }

   return erg;
}




static int gives_rise_to_torsionfree_space_group(
         bravais_TYP *G,
         matrix_TYP *D,
         matrix_TYP *cocycle,
         matrix_TYP *x)
{

    bravais_TYP *R;      /* holds the space group */

    int i,
        j,
        k,
        *res;

    matrix_TYP *C;


    R = init_bravais(G->dim+1);
    C = convert_to_cozycle(x,cocycle,D);

    R->gen = (matrix_TYP **) malloc(G->gen_no * sizeof(matrix_TYP *));
    R->gen_no = G->gen_no;
    for (i=0;i<G->gen_no;i++){
       R->gen[i] = copy_mat(G->gen[i]);
       real_mat(R->gen[i],G->dim+1,G->dim+1);
       R->gen[i]->array.SZ[G->dim][G->dim] = 1;
       iscal_mul(R->gen[i],C->kgv);
       R->gen[i]->kgv = C->kgv;
    }

    /* stick the cocycle in the last column of the matrices generating R */
    k = 0;
    for (i=0;i<R->gen_no;i++){
       for (j=0;j<G->dim;j++){
          R->gen[i]->array.SZ[j][G->dim] = C->array.SZ[k][0];
          k++;
       }
       Check_mat(R->gen[i]);
    }

    res = torsionfree(R,&i,&j);
    i = res[0] ; free(res);

    free_bravais(R);
    free_mat(C);

    return i;
}




/**************************************************************************
@
@ -------------------------------------------------------------------------
@
@ matrix_TYP **extensions(matrix_TYP *cozycle,
@                         matrix_TYP *D,
@                         matrix_TYP *R,
@                         bravais_TYP *G,
@                         int **lengths,
@                         MP_INT **names,
@                         int *number_of_orbits,
@                         int option)
@
@ Returns the cozycles which generate the isomorphims classes of
@ extensions of G by the natural ZG-module. The split extension is
@ represented by an all zero matrix.
@
@ The arguments are:
@   matrix_TYP *cozykle:   1st matrix returned by cohomolgy for G.
@   matrix_TYP *D:         2nd matrix returned by cohomolgy for G.
@   matrix_TYP *R:         3rd matrix returned by cohomolgy for G.
@   bravais_TYP *G:        the group in question.
@   int **lengths:         length[0] returns a pointer to the lengths
@                          of the orbits respectively
@   MP_INT **names:        names[0] returns a pointer to the names of
@                          the cozycles as they would appear in a call
@                          of identify(.....).
@   int *number_of_orbits: the number of orbits the normalizer induces
@                          on the cohomology group.
@   MP_INT coho_size:      order of the cohomology group
@   int option:            controls the behaviour of the function:
@                          option & 1: construct only torsion free extensions
@   int *list_of_names:    if list_of_names != NULL save of each cocycle the number of the
@                          smallest represenatative in it's orbit there
@ -------------------------------------------------------------------------
@
***********************************************************************/
matrix_TYP **extensions_o(matrix_TYP *cozycle,
                          matrix_TYP *D,
                          matrix_TYP *R,
                          bravais_TYP *G,
                          int **lengths,
                          MP_INT **names,
                          int *number_of_orbits,
                          int ****WORDS,
                          int **NUMBER_OF_WORDS,
                          matrix_TYP ***N,
                          MP_INT coho_size,
                          int option,
                          int *list_of_names)
{
  int i, j, k, anz, *wort, act_val_int, coho_size_int,
      Nanz = G->normal_no + G->cen_no,
      orbit_length, **WWW;

  char *tested;

  MP_INT act_val,
         new_val,
         got_size;

  matrix_TYP *x,
             *y,
            **erg;  /* stores the representatives of the orbits=isomorphism
                       classes of extensions as rational matrices */

  if (INFO_LEVEL & 7){
     fprintf(stderr,"entered extensions\n");
  }

  number_of_orbits[0] = 0;  /* should be changed to a calculation via burnside's
                               lemma */

  mpz_init_set_si(&act_val,0);
  mpz_init_set_si(&got_size,0);
  mpz_init(&new_val);

  erg = (matrix_TYP **) malloc(sizeof(matrix_TYP *));
  lengths[0] = (int *) malloc(sizeof(int));
  names[0] = (MP_INT *) malloc(sizeof(MP_INT));
  N[0] = (matrix_TYP **)calloc(Nanz, sizeof(matrix_TYP *));

  if (list_of_names != NULL){
     coho_size_int = mpz_get_ui(&coho_size);
  }

  if (mpz_get_ui(&coho_size) < TWOTO21){
     tested = (char *)calloc(mpz_get_ui(&coho_size), sizeof(char));
  }
  else{
     tested = NULL;
     if (list_of_names != NULL){
        fprintf(stderr, "The cohomology group is too big!\n");
        fprintf(stderr, "If you are really interested in the graph of group - subgroup relations\n");
        fprintf(stderr, "please start the program again with the option -l!\n");
        exit(3);
     }
  }

  for (i=0;i<Nanz;i++){
     if (i<G->cen_no){
        N[0][i] = normalop(cozycle,D,R,G,G->cen[i],i == (Nanz-1));
     }
     else{
        N[0][i] = normalop(cozycle,D,R,G,G->normal[i-G->cen_no],i == (Nanz-1));
     }
     if (INFO_LEVEL & 16){
        put_mat(N[0][i],NULL,"N[i]",2);
     }
  }

  if (is_option('v')){
     fprintf(stderr,"The 1-cohomology group has ");
     mpz_out_str(stderr,10,&coho_size);
  }

  while((mpz_cmp(&got_size, &coho_size)<0) && (mpz_cmp(&act_val, &coho_size)<0)){

     if (INFO_LEVEL == 17 ||
         is_option('v')){
        fprintf(stderr,"Got %d extensions so far\n",*number_of_orbits);
     }

     if (tested == NULL || !tested[mpz_get_ui(&act_val)]){

        if (tested != NULL) tested[mpz_get_ui(&act_val)] = TRUE;

        WWW = (int **)calloc(MIN_SPEICHER, sizeof(int *));
        anz = 0;
        x = reverse_valuation(&act_val,D);
        y = static_orbit_rep(x,N[0],Nanz,D,1,tested,&new_val,&orbit_length,
                             &wort, 1, &WWW, &anz);
        free(wort);

        if ((INFO_LEVEL & 4) || list_of_names != NULL)
           act_val_int = mpz_get_ui(&act_val);

        if (INFO_LEVEL & 4){
           printf("act_val %d\n", act_val_int);
           put_mat(x,NULL,"x",2);
           put_mat(y,NULL,"y",2);
        }

        if (list_of_names != NULL){
           for (k = 1; k < coho_size_int; k++){
              if (list_of_names[k] == 0 && tested[k] == TRUE)
                 list_of_names[k] = act_val_int;
           }
        }

        if (mpz_cmp(&new_val,&act_val) == 0 &&
           (!(option & 1) ||
              gives_rise_to_torsionfree_space_group( G, D, cozycle, x))){
          if (number_of_orbits[0] % MIN_SPEICHER == 0){
             erg = (matrix_TYP **) realloc(erg,
                     (number_of_orbits[0] + MIN_SPEICHER)*sizeof(matrix_TYP*));
             lengths[0] = (int *) realloc(lengths[0],
                     (number_of_orbits[0] + MIN_SPEICHER)*sizeof(int));
             names[0] = (MP_INT *) realloc(names[0],
                     (number_of_orbits[0] + MIN_SPEICHER)*sizeof(MP_INT));
             WORDS[0] = (int ***)realloc(WORDS[0],
                     (number_of_orbits[0] + MIN_SPEICHER)*sizeof(int **));
             NUMBER_OF_WORDS[0] = (int *)realloc(NUMBER_OF_WORDS[0],
                     (number_of_orbits[0] + MIN_SPEICHER)*sizeof(int));

          }
          erg[number_of_orbits[0]] = convert_to_cozycle(x,cozycle,D);
          lengths[0][number_of_orbits[0]] = orbit_length;
          mpz_init_set(names[0]+number_of_orbits[0],&act_val);
          mpz_add_ui(&got_size,&got_size,(unsigned long) orbit_length);
          WORDS[0][number_of_orbits[0]] = (int **)calloc(anz, sizeof(int *));
          for (j = 0; j < anz; j++){
             WORDS[0][number_of_orbits[0]][j] = WWW[j];
          }
          NUMBER_OF_WORDS[0][number_of_orbits[0]] = anz;
          free(WWW);
          number_of_orbits[0]++;
        }
        else{
           for (j = 0; j < anz; j++)
              free(WWW[j]);
           free(WWW);
        }

        /* clean up and increase act_val */
        free_mat(y);
        free_mat(x);
     }
     mpz_add_ui(&act_val,&act_val,1);
  }

  if (tested != NULL ) free(tested);
  mpz_clear(&act_val);
  mpz_clear(&new_val);
  mpz_clear(&got_size);

  return erg;
}  /* extensions (.... ) */




/***************************************************************************/
void my_translation(matrix_TYP *TR,
                    matrix_TYP *preimage,
                    matrix_TYP *coz,
                    bravais_TYP *G)
{

   int i,
       j,
       k,
       **WORDS;

   matrix_TYP *tmp,
              *tmp2,
             **MAP;

   rational eins, minuseins;


   eins.z = eins.n = minuseins.n = 1; minuseins.z = -1;
   WORDS = (int **)calloc(G->gen_no, sizeof(int *));

   tmp2 = long_mat_inv(TR);
   MAP = (matrix_TYP **) malloc(G->gen_no * sizeof(matrix_TYP *));
   for (j=0;j<G->gen_no;j++){
      MAP[j] = copy_mat(G->gen[j]);
      real_mat(MAP[j],MAP[j]->rows+1,MAP[j]->cols);
      real_mat(MAP[j],MAP[j]->rows,MAP[j]->cols+1);
      MAP[j]->array.SZ[G->dim][G->dim] = 1;
      for (k=0;k<G->dim;k++){
         MAP[j]->array.SZ[k][G->dim]=preimage->array.SZ[j*G->dim+k][0];
      }
      tmp = mat_kon(TR,MAP[j],tmp2);
      free_mat(MAP[j]);
      MAP[j] = tmp;
   }
   tmp = reget_gen(MAP,G->gen_no,G,WORDS,TRUE);
   tmp->kgv = preimage->kgv;

   for (i=0;i<G->gen_no;i++){
      free_mat(MAP[i]);
      free(WORDS[i]);
   }
   free(MAP);
   free(WORDS);

   tmp = mat_addeq(tmp, coz, eins, minuseins);

   coboundary(G,tmp,TR);
   free_mat(tmp2);
   free_mat(tmp);

   return;

}


