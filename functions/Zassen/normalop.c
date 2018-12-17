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
#include <datei.h>

#define TWOTO21 2097152

extern int INFO_LEVEL;

static void mod(int *a,int b)
{
   a[0] = a[0] % b;
   if (a[0]<0) a[0] += b;
}

void valuation(matrix_TYP *x,matrix_TYP *D,MP_INT *val)
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
matrix_TYP *orbit_rep(matrix_TYP *x,
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
   valuation(x,D,&valuations[0]);

   /* anz[0] is the length of the orbit so far */
   anz[0] = 1;
   erg = NULL;

   for (i=0;i<anz[0] && erg == NULL;i++){
      for (j=0;j<nanz && erg == NULL;j++){
         tmp = mat_mul(N[j],orbit[i]);

         /* normalize tmp modulo the diagonal entries of D */
         for (k=first;k<last;k++)
           mod(tmp->array.SZ[k-first],D->array.SZ[k][k]);

         valuation(tmp,D,&new_val);

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



/**************************************************************************
@
@ -------------------------------------------------------------------------
@
@ matrix_TYP *normalop(matrix_TYP *cozycle,matrix_TYP *D,matrix_TYP *R,
@                      bravais_TYP *G,matrix_TYP *N)
@
@ Calculates the action of the matrix N on the cohomology group H^1(G,*)
@ described by the matrices cozycle, D, R, which are in turn output of
@ cohomology.
@
@ CAUTION: The matrix returned describes the action on rows, and it is
@          not checked whether N is realy an element of the normalizer of G!
@
@ matrix_TYP *cozycle: the first return value of cohomology for G.
@ matrix_TYP *D:       the second return value of cohomology for G.
@ matrix_TYP *R:       the third return value of cohomology for G.
@ bravais_TYP *G:      the group. only the field generator is realy needed.
@ matrix_TYP  *N:      the matrix in question.
@ int opt:
@
@ No global variables nor the arguments are changed (hopefully).
@
@ -------------------------------------------------------------------------
@
***************************************************************************/
matrix_TYP *normalop(matrix_TYP *cozycle,
                     matrix_TYP *D,
                     matrix_TYP *R,
                     bravais_TYP *G,
                     matrix_TYP *N,
                     int opt)

{
  int i,
      j,
      k,
      denominator,
      first,      /* first index such that D[i][i] != 1 */
      last;       /* last index with D[i][i] !=0 */

   int **words;  /* the words which give the original generators
                    will be stored in here by reget_gen */

   matrix_TYP *tmp,
              *tmp2,
              *sols,
              *new_coz, /* the mapped cozycle */
              *Ninv,    /* hold's the inverse of N */
             **map,     /* hold's the map of a pair generator/cozycle under N */
              *erg;     /* describes the result */

   static matrix_TYP *GLS;

   Ninv = mat_inv(N);
   erg = init_mat(cozycle->cols,cozycle->cols,"");
   map = (matrix_TYP **) malloc(G->gen_no * sizeof(matrix_TYP *));
   for (i=0;i<G->gen_no;i++){
      map[i] = init_mat(G->dim+1,G->dim+1,"1");
   }

   /* get memory for the words */
   words = (int **) malloc(G->gen_no * sizeof(int *));

   /* set first and last */
   for (first = 0;first<D->cols && D->array.SZ[first][first] == 1;first++);
   for (last = 0;last<D->cols && D->array.SZ[last][last] != 0;last++);

   if (GLS == NULL){
      GLS = long_mat_inv(R);
   }

   if (INFO_LEVEL & 4){
      fprintf(stderr,"entered normalop\n");
      put_mat(N,NULL,"N",2);
      put_mat(Ninv,NULL,"Ninv",2);
      put_mat(GLS,NULL,"GLS",2);
   }

   /* we take a generating set of the cozycles and let this
      normalizer element operate */
   for (i=first;i<last;i++){
      /* map the i-th cozycle */
      for (j=0;j<G->gen_no;j++){
         /* map the jth generator */
         tmp = mat_kon(N,G->gen[j],Ninv);

         /* store it in the right place */
         for (k=0;k<tmp->rows;k++){
            memcpy(map[j]->array.SZ[k],tmp->array.SZ[k],
                   sizeof(int) * tmp->cols);
         }

         free_mat(tmp);

         /* now map the cozycle */
         tmp = init_mat(G->dim,1,"");
         for (k=0;k<tmp->rows;k++)
            tmp->array.SZ[k][0] = cozycle->array.SZ[j*G->dim + k][i-first];

         tmp2 = mat_mul(N,tmp);
         free_mat(tmp);

         for (k=0;k<tmp->rows;k++)
            map[j]->array.SZ[k][G->dim] = tmp2->array.SZ[k][0];
         free_mat(tmp2);

         if (INFO_LEVEL & 16){
            Check_mat(map[j]);
            printf("map of the %d-th gen/coz\n",j);
            put_mat(map[j],NULL,"map[j]",2);
         }
      }

      /* reget the original generators of the group,
         bare in mind that we got words for i>first already */
      new_coz = reget_gen(map,G->gen_no,G,words,i==first);

      /* now we are able to express this new cozycle in terms of the "basis" */
      /* solve the resulting linear system of equations */
      sols = mat_mul(GLS,new_coz);

      if (INFO_LEVEL & 16){
         put_mat(sols,NULL,"sols",2);
         put_mat(new_coz,NULL,"new_coz",2);
      }

      /* copy the solution to the COLUMNS of erg, so erg will
         act on COLUMNS as well */
      for (k=0;k<erg->rows;k++){
         j = sols->array.SZ[k+first][0] * D->array.SZ[k+first][k+first];
         denominator = D->array.SZ[i][i];
         if ((j % denominator)!=0){
            fprintf(stderr,"error in normalop\n");
            exit(3);
         }
         else{
            erg->array.SZ[k][i-first] = (j / denominator) %
                                         D->array.SZ[k+first][k+first];
            if (erg->array.SZ[k][i-first] < 0)
               erg->array.SZ[k][i-first] += D->array.SZ[k+first][k+first];
         }
      }

      free_mat(sols);
      free_mat(new_coz);
      
   }

   /* cleaning up */
   free_mat(Ninv);
   for (i=0;i<G->gen_no;i++){
      free_mat(map[i]);
      free(words[i]);
   }
   free(map);
   free(words);

   if (opt){
      free_mat(GLS);
      GLS = NULL;
   }

   return erg;
}


void translation(matrix_TYP *TR,
                 matrix_TYP *rep,
                 matrix_TYP *ext,
                 matrix_TYP *cocycle,
                 matrix_TYP *D,
                 bravais_TYP *G)
{

   int i,
       j,
       k,
       first,
     **WORDS;

   matrix_TYP *tmp,
              *tmp2,
             **MAP;

   /* set first */
   for (first = 0;first<D->cols && D->array.SZ[first][first] == 1;first++);

   WORDS = (int **) calloc(G->gen_no , sizeof(int*));

   tmp2 = long_mat_inv(TR);
   MAP = (matrix_TYP **) malloc(G->gen_no * sizeof(matrix_TYP *));
   for (j=0;j<G->gen_no;j++){
      MAP[j] = copy_mat(G->gen[j]);
      real_mat(MAP[j],MAP[j]->rows+1,MAP[j]->cols);
      real_mat(MAP[j],MAP[j]->rows,MAP[j]->cols+1);
      MAP[j]->array.SZ[G->dim][G->dim] = 1;
      for (k=0;k<G->dim;k++){
         MAP[j]->array.SZ[k][G->dim]=ext->array.SZ[j*G->dim+k][0];
      }
      tmp = mat_kon(TR,MAP[j],tmp2);
      free_mat(MAP[j]);
      MAP[j] = tmp;
   }
   tmp = reget_gen(MAP,G->gen_no,G,WORDS,TRUE);
   tmp->kgv = ext->kgv;

   for (i=0;i<G->gen_no;i++){
      free_mat(MAP[i]);
      free(WORDS[i]);
   }
   free(MAP);
   free(WORDS);

   for (i=0;i<rep->rows;i++){
      k = rep->array.SZ[i][0] * tmp->kgv / D->array.SZ[first+i][first+i];
      if (k != 0){
         for (j=0;j<tmp->rows;j++){
            tmp->array.SZ[j][0] -= k * cocycle->array.SZ[j][i];
         }
      }
   }

   coboundary(G,tmp,TR);
   free_mat(tmp2);
   free_mat(tmp);

   return;

}

/**************************************************************************
@
@--------------------------------------------------------------------------
@
@ matrix_TYP **identify(matrix_TYP *cozycle,
@                       matrix_TYP *D,
@                       matrix_TYP *R,
@                       bravais_TYP *G,
@                       matrix_TYP **extension,
@                       MP_INT *a,
@                       int number,
@                       int transform_flag,
@                       int ***WORDS,
@                       int *NUMBER_OF_WORDS)
@
@ Identifies the space groups described by the cozycles in extension,
@ i.e. gives them different number iff they are not isomorphic.
@ The return value depends on the value of transform_flag. If
@ (transform_flag == TRUE) it will return a list of "number" matrices in
@ the normalizer of the point group which transforms each given extension
@ into it's least representative. This is usefull to calculate the
@ isomorphism between two extensions which get the same name.
@ NOTE: the name given to an extension is 0 (interpreted as MP_INT) iff
@       the extension splits.
@
@   matrix_TYP *cozykle:   1st matrix returned by cohomolgy for G.
@   matrix_TYP *D:         2nd matrix returned by cohomolgy for G.
@   matrix_TYP *R:         3rd matrix returned by cohomolgy for G.
@   bravais_TYP *G:        the group in question. Needed are the
@                          fields G->gen, G->gen_no, G->cen,G->cen_no,
@                          G->normal,G->normal_no.
@                          Important for the correctness of the result
@                          is that the matrices in G->cen, G->normal
@                          generate N_GL_n(Z) (G) /G.
@  matrix_TYP **extension: matrices discribing a cozycles for the group G.
@                          So each matrix in the list discribes a spacegroup.
@  MP_int *a:              The names for the groups will be returned
@                          via this field. The space has to be allocated
@                          via malloc and to be mpz_init-ed for each
@                          entry before calling the function.
@  int number:             The number of matrices in extension.
@  int transform_flag:     if (transform_flag) the function will return
@                          a pointer containing "number" matrices which
@                          transform each extension into it's least
@                          representative. Otherwise it will return NULL.
@  int ***WORDS:           words in the generators of the normalizer
@                          which generate the stabilizer og this cozykle
@                          as a subgroup. If == NULL, this is not calculated.
@                          Otherwise it is an address of a dynamic pointer
@                          to at least MIN_SPEICHER entries of type *int.
@  int *NUMBER_OF_WORDS:   number of words WORDS.
@--------------------------------------------------------------------------
@
***************************************************************************/
matrix_TYP **identify(matrix_TYP *cozycle,
                      matrix_TYP *D,
                      matrix_TYP *R,
                      bravais_TYP *G,
                      matrix_TYP **extension,
                      MP_INT *a,
                      int number,
                      int transform_flag,
                      int ***WORDS,
                      int *NUMBER_OF_WORDS)
{

  int i,
      j,
      k,
      denominator,
      first,
     *word,                        /* the function orbit_rep will return
                                      a word representing the least extension
                                      by this pointer */
      Nanz= G->cen_no+G->normal_no;

  matrix_TYP **N,
             **TR,          /* holds the result, ie. the transformation
                               matrices */
              *GLS,
              *tmp,
              *coz,
              *rep;

  /* set first */
  for (first = 0;first<D->cols && D->array.SZ[first][first] == 1;first++);

  if (transform_flag)
      TR = (matrix_TYP **) malloc(number * sizeof(matrix_TYP *));
  else
      TR = NULL;

  GLS = mat_inv(R);
  N = (matrix_TYP **) malloc(Nanz * sizeof(matrix_TYP*));
  for (i=0;i<Nanz;i++){
     if (i<G->cen_no){
        N[i] = normalop(cozycle,D,R,G,G->cen[i],i == (Nanz-1));
     }
     else{
        N[i] = normalop(cozycle,D,R,G,G->normal[i-G->cen_no],i == (Nanz-1));
     }
     if (INFO_LEVEL & 16){
        put_mat(N[i],NULL,"N[i]",2);
     }
  }

  coz = init_mat(N[0]->rows,1,"");

  for (i=0;i<number;i++){

     /* firstly calculate a representative, and translate it into
        the representation of C_d1 X ..... X C_dn */
     tmp = mat_mul(GLS,extension[i]);

     /* and now calculate the smallest representative in the orbit under
        the action of the normalizer, and it's valuation */
     for (k=0;k<coz->rows;k++){
        j = tmp->array.SZ[k+first][0] * D->array.SZ[k+first][k+first];
        denominator = tmp->kgv;
        if ((j % denominator)!=0){
           fprintf(stderr,"error in identify: are you sure this is a cozycle?\n");
           fprintf(stderr,"If so, please report to the authors: carat@momo.math.rwth-aachen.de\n");
           exit(3);
        }
        else{
           coz->array.SZ[k][0] = (j / denominator) %
                                        D->array.SZ[k+first][k+first];
           while (coz->array.SZ[k][0] < 0)
              coz->array.SZ[k][0] += D->array.SZ[k+first][k+first];
        }
     }

     rep = orbit_rep(coz,N,Nanz,D,0,NULL,a+i,&j,&word,
                       transform_flag,WORDS,NUMBER_OF_WORDS);

     if (transform_flag){
        TR[i] = init_mat(G->dim,G->dim,"1");
        for (j=1;j<=word[0];j++){
          if (word[j] < G->cen_no){
             mat_muleq(TR[i],G->cen[word[j]]);
          }
          else{
             mat_muleq(TR[i],G->normal[word[j]-G->cen_no]);
          }
        }
        free(word);

        /* now calculate the coboundary to get it physically nailed */
        real_mat(TR[i],TR[i]->rows,TR[i]->cols+1);
        real_mat(TR[i],TR[i]->rows+1,TR[i]->cols);
        TR[i]->array.SZ[G->dim][G->dim] = 1;

        if (transform_flag == 3){
           translation(TR[i],rep,extension[i],cozycle,D,G);
        }
     }

     /* we are done, so clean up */
     free_mat(tmp);
     free_mat(rep);
     
  }

  /* clean up N, coz and GLS */
  for (i=0;i<Nanz;i++) free_mat(N[i]);
  free(N);
  free_mat(coz);
  free_mat(GLS);

  return TR;
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
@    int option          : controls the behaviour of the function:
@                          option & 1: construct only torsion free extensions
@ -------------------------------------------------------------------------
@
***************************************************************************/
matrix_TYP **extensions(matrix_TYP *cozycle,
                        matrix_TYP *D,
                        matrix_TYP *R,
                        bravais_TYP *G,
                        int **lengths,
                        MP_INT **names,
                        int *number_of_orbits,
                        int option )
{
  int i,
      Nanz = G->normal_no + G->cen_no,
      orbit_length;

  char *tested;

  MP_INT act_val,
         new_val,
         got_size,
         coho_size;   /* holds the size of the cohomology group */

  matrix_TYP *x,
             *y,
            **N,
            **erg;  /* stores the representatives of the orbits=isomorphism
                       classes of extensions as rational matrices */

  if (INFO_LEVEL & 7){
     fprintf(stderr,"entered extensions\n");
  }

  number_of_orbits[0] = 0;  /* should be changed to a calculation via burnside's
                               lemma */

  mpz_init_set_si(&act_val,0);
  mpz_init_set_si(&coho_size,1);
  mpz_init_set_si(&got_size,0);
  mpz_init(&new_val);
  /* the order of the cohomology group in the product of those elementary
     divisors which are not equal 0 */
  for (i=0;i<D->cols && D->array.SZ[i][i] != 0;i++)
           mpz_mul_ui(&coho_size,&coho_size,(unsigned long) D->array.SZ[i][i]);

  erg = (matrix_TYP **) malloc(sizeof(matrix_TYP *));
  lengths[0] = (int *) malloc(sizeof(int));
  names[0] = (MP_INT *) malloc(sizeof(MP_INT));
  N = (matrix_TYP **) malloc(Nanz * sizeof(matrix_TYP *));

  if (mpz_get_ui(&coho_size) < TWOTO21)
     tested = (char *) calloc(mpz_get_ui(&coho_size), sizeof(char));
  else
     tested = NULL;

  for (i=0;i<Nanz;i++){
     if (i<G->cen_no){
        N[i] = normalop(cozycle,D,R,G,G->cen[i],i == (Nanz-1));
     }
     else{
        N[i] = normalop(cozycle,D,R,G,G->normal[i-G->cen_no],i == (Nanz-1));
     }
     if (INFO_LEVEL & 16){
        put_mat(N[i],NULL,"N[i]",2);
     }
  }

  if (is_option('v')){
     fprintf(stderr,"The 1-cohomology group has ");
     mpz_out_str(stderr,10,&coho_size);
  }

  while((mpz_cmp(&got_size,&coho_size)<0) && (mpz_cmp(&act_val,&coho_size)<0)){

     if (INFO_LEVEL == 17 ||
         is_option('v')){
        fprintf(stderr,"Got %d extensions so far\n",*number_of_orbits);
     }

     if (tested == NULL || !tested[mpz_get_ui(&act_val)]){

        if (tested != NULL) tested[mpz_get_ui(&act_val)] = TRUE;

        x = reverse_valuation(&act_val,D);
        y = orbit_rep(x,N,Nanz,D,1,tested,&new_val,&orbit_length,
                             NULL,FALSE,NULL,NULL);

        if (INFO_LEVEL & 4){
           printf("act_val %d\n",(int ) mpz_get_ui(&act_val));
           put_mat(x,NULL,"x",2);
           put_mat(y,NULL,"y",2);
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
          }
          erg[number_of_orbits[0]] = convert_to_cozycle(x,cozycle,D);
          lengths[0][number_of_orbits[0]] = orbit_length;
          mpz_init_set(names[0]+number_of_orbits[0],&act_val);
          number_of_orbits[0]++;
          mpz_add_ui(&got_size,&got_size,(unsigned long) orbit_length);
        }

        /* clean up and increase act_val */
        free_mat(y);
        free_mat(x);
     }
     mpz_add_ui(&act_val,&act_val,1);
  }

  for (i=0;i<Nanz;i++) free_mat(N[i]);
  free(N);
  if (tested != NULL ) free(tested);
  mpz_clear(&act_val);
  mpz_clear(&new_val);
  mpz_clear(&got_size);
  mpz_clear(&coho_size);

  return erg;
}  /* extensions (.... ) */

static matrix_TYP **calc_elements(matrix_TYP **N,matrix_TYP **Ninv,int Nanz,
              matrix_TYP *D,unsigned long *no_of_elements)
{
  int i,
      j,
      k,
      l,
      pos,        /* position returned by hash_mat */
      dim,        /* the rank of the Z-module the N[i] are acting on */
      first,      /* first index such that D[i][i] != 1 */
      last;       /* last index with D[i][i] !=0 */

   matrix_TYP  *tmp,
              **erg;      /* holds all elements of the group */

   struct tree *hash;

   /* set first and last */
   for (first = 0;first<D->cols && D->array.SZ[first][first] == 1;first++);
   for (last = 0;last<D->cols && D->array.SZ[last][last] != 0;last++);
   dim = last-first;

   erg = (matrix_TYP **) malloc(MIN_SPEICHER * sizeof(matrix_TYP *));
   erg[0] = init_mat(dim,dim,"1");
   *no_of_elements = 1;
   hash = (struct tree *) malloc(1*sizeof(struct tree));
   hash->no = 0;
   hash->left = hash->right = NULL;

   /* we won't need inverses for the orbit calculation, everything is finite */
   for (i=0;i<*no_of_elements;i++){
      for (j=0;j<Nanz;j++){
         tmp = mat_mul(erg[i],N[j]);

         /* normalize tmp */
         for (k=0;k<dim;k++)
            for (l=0;l<dim;l++)
               mod(tmp->array.SZ[k]+l,D->array.SZ[k+first][k+first]);

         pos = hash_mat(hash,erg,tmp,*no_of_elements);

         if (pos != -1){
            if (Ninv!=NULL && pos==0) Ninv[j] = erg[i];
            free_mat(tmp);
         }
         else{
            if (*no_of_elements % MIN_SPEICHER == 0)
               erg = (matrix_TYP **) realloc(erg,
                      (*no_of_elements + MIN_SPEICHER)*sizeof(matrix_TYP *));

            erg[*no_of_elements] = tmp;
            no_of_elements[0]++;
         }
      }
   }

   free_tree(hash);

   return erg;
} /* calc_elements(....) */


/****************************************************************************
@
@----------------------------------------------------------------------------
@
@ static void no_of_fixpoints(MP_INT *res,matrix_TYP *A,matrix_TYP *D)
@
@
@ CAUTION: A->array.SZ, A->cols WILL BE CHANGED!!!!
@----------------------------------------------------------------------------
@
*****************************************************************************/
static void no_of_fixpoints(MP_INT *res,matrix_TYP *A,matrix_TYP *D)
{
   int i,
       first;

   matrix_TYP *tmp;

   /* set first */
   for (first = 0;first<D->cols && D->array.SZ[first][first] == 1;first++);

  /* build a matrix (A-I | D') in A (D' = those parts of D which have non
     trivial elementary divisors) */
  real_mat(A,A->rows,A->cols+A->rows);
  for (i=0;i<A->rows;i++){
     A->array.SZ[i][i]--;
     A->array.SZ[i][i+A->rows] = D->array.SZ[i+first][i+first];
  }

  tmp = long_elt_mat(NULL,A,NULL);

  mpz_set_ui(res,1);
  for (i=0;i<tmp->rows;i++) mpz_mul_ui(res,res,tmp->array.SZ[i][i]);

  free_mat(tmp);
}

static int deal_with_small_cohomology_group(matrix_TYP *cozycle,
                                            matrix_TYP *D,
                                            matrix_TYP *R,
                                            bravais_TYP *G,
                                            MP_INT *erg)
{

  int i,
      anz,
      order_of_H=1,
     *len;

  MP_INT *names;

  matrix_TYP **Y;

  if (G->dim < 5){
     return FALSE;
  }

  /* if the cohomology is small, it might be better to calculate all
     extensions then to calculate all elements of the acting group
     example: H^1 = C_6^3, G_acting = GL_3(Z/6Z) */
  for (i=0;i < D->cols &&
           i < D->rows &&
           order_of_H < 10000 &&
           D->array.SZ[i][i] != 0;i++,order_of_H *= D->array.SZ[i][i]);
  if (order_of_H < 10000){
     Y = extensions(cozycle,D,R,G,&len,&names,&anz,0);

     for (i=0;i<anz;i++){
        free_mat(Y[i]);
        mpz_clear(names+i);
     }
     free(names);
     free(Y);
     free(len);
     mpz_set_si(erg,anz);
     return TRUE;
  }

  return FALSE;

}

/************************************************************************
@
@------------------------------------------------------------------------
@
@ void no_of_extensions(matrix_TYP *cozycle,matrix_TYP *D,
@                        matrix_TYP *R,bravais_TYP *G,MP_INT *erg)
@
@------------------------------------------------------------------------
@
*************************************************************************/
void no_of_extensions(matrix_TYP *cozycle,matrix_TYP *D,
                        matrix_TYP *R,bravais_TYP *G,MP_INT *erg)
{

  int Nanz=G->cen_no+G->normal_no,
      i;

  unsigned long no_of_elements;   /* the order of the image of the
                                     normalizer under the map to the
                                     automorphism group of H^1(...) */

  MP_INT sum,       /* holds the sum over the image of N_GL(Z) in
                       Aut(H^1(...)) of the number fixpoint */
         tmp2,
         tmp3;

  matrix_TYP **N,
             **elements,
              *id;

  if (deal_with_small_cohomology_group(cozycle, D, R, G, erg)){
     return;
  }

  mpz_init(&sum);
  mpz_set_si(&sum,0);
  mpz_init(&tmp2);
  mpz_set_si(&tmp2,0);
  mpz_init(&tmp3);
  mpz_set_si(&tmp3,0);

  N = (matrix_TYP **) malloc(Nanz * sizeof(matrix_TYP *));
  for (i=0;i<Nanz;i++){
     if (i<G->cen_no){
        N[i] = normalop(cozycle,D,R,G,G->cen[i],i == (Nanz-1));
     }
     else{
        N[i] = normalop(cozycle,D,R,G,G->normal[i-G->cen_no],i == (Nanz-1));
     }
     if (INFO_LEVEL & 16){
        put_mat(N[i],NULL,"N[i]",2);
     }
  }

  /* remove those generators which are == identity for speeding up things */
  if (Nanz > 0){
    id = init_mat(N[0]->rows,N[0]->rows,"1");
    for (i=0;i<Nanz;i++){
       if (mat_comp(id,N[i]) == 0){
          free_mat(N[i]);
          Nanz--;
          N[i] = N[Nanz];
       }
    }
  }

  if (is_option('v')){
     fprintf(stderr,"calculated the action of the normalizer on the\n");
     fprintf(stderr,"cohomology group\n");
  }

  elements = calc_elements(N,NULL,Nanz,D,&no_of_elements);

  if (is_option('v')){
     fprintf(stderr,"Now I got all elements in the image.\n");
  }

  for (i=0;i<no_of_elements;i++){
    no_of_fixpoints(&tmp2,elements[i],D);

    if (INFO_LEVEL & 4){
       mpz_out_str(stdout,10,&tmp2);
       printf("\n");
    }

    mpz_add(&tmp3,&sum,&tmp2);
    mpz_set(&sum,&tmp3);
    /* we won't need this element anymore */
    free_mat(elements[i]);

    if (is_option('v') && (i % 10 == 0)){
       fprintf(stderr,"Done approximately %d percent of the calculation\n",
                       (int ) ((i*100)/no_of_elements));
    }

  }

  /* divide sum by the image order, and check whether this is possible
     at all */
  mpz_divmod_ui(&tmp3,&tmp2,&sum,no_of_elements);
  mpz_set(erg,&tmp3);

  if (INFO_LEVEL & 4){
     mpz_out_str(stdout,10,&sum);
     printf("\n");
     mpz_out_str(stdout,10,erg);
     printf("\n");
     mpz_out_str(stdout,10,&tmp2);
     printf("\n");
  }

  if (mpz_cmp_si(&tmp2,0) != 0){
      printf("error in no_of_extensions\n");
      exit(3);
  }

  for (i=0;i<Nanz;i++) free_mat(N[i]);
  free(N);
  free(elements);
  mpz_clear(&sum);
  mpz_clear(&tmp2);
  mpz_clear(&tmp3);

  return;
} /* no_of_extensions(....) */
