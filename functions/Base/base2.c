/***************************************************************************
@
@ --------------------------------------------------------------------------
@
@ FILE: base2.c
@
@ some algorithms connected with the schreier-sims algorithm are collected.
@ in contrast to the functions in base.c these functions operate mod 2,
@ so these functions are able to calculate the kernel of the minkowski
@ homomorphism as well.
@
@ --------------------------------------------------------------------------
@
****************************************************************************/

#include <typedef.h>
#include <longtools.h>
#include <orbit.h>
#include <getput.h>
#include <matrix.h>
#include <sort.h>
#include <base.h>


/***************************************************************************
@
@ --------------------------------------------------------------------------
@
@  int einordnen_2(bahn** erg,matrix_TYP *h, matrix_TYP *new_vec,int l,
@                       int flag, matrix_TYP ***K,int *anz)
@
@  bahn **erg          :  the set of base/strong generators (mod 2) for a
@                         (possibly infinite) group G. this variable will
@                         be changed!
@  matrix_TYP *h       :
@  matrix_TYP *new_vec :  this value is not asked for in the case
@                         described here, so any variable with the right
@                         type will do.
@  int l               :  -1          ALLWAYS
@  int flag            :  TRUE        ALLWAYS!!!!!
@  matrix_TYP ***K     :  this list of matrices might already contain
@                         generators for the kernel of the minkowski
@                         homomorphism. All newly found generators
@                         will be stucked in K[0][anz[0]],K[0][anz[0]+1]..
@                         and anz[0] will be changed accordingly.
@  int *anz[0]         :  see above
@
@ Inserts a new generator h into the list of strong generators (mod 2) for
@ any (possible infinite) group G.
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
int einordnen_2(bahn** erg,matrix_TYP *h, matrix_TYP *new_vec,int l,
                      int flag,matrix_TYP ***K,int *anz)
{
   int tmp = erg[0]->representatives[0]->cols,
       tmp2,
       i,
       finite_flag = TRUE;

   matrix_TYP  *hinv,
               *nv,
               *nh;

   if (INFO_LEVEL & 2){
      fprintf(stderr,"entering einordnen_2 with l = %d\n",l);
   }

   if (l>=tmp){
      /* this is the stage to memorize the kernel of
         the minkowski homomorphism */
      for (l=0;l<anz[0] && mat_comp(K[0][l],h) != 0;l++);

      if (anz[0] == 0 || l == anz[0]){
        /* get more memory if needed */
        if (anz[0]>0 && (anz[0] % MIN_SPEICHER) == 0){
           K[0] = (matrix_TYP **) realloc(K[0],
                              (anz[0]+MIN_SPEICHER)*sizeof(matrix_TYP *));
        }
        K[0][anz[0]] = h;
        anz[0]++;

        /* check if the group stay's finite */
        /* the new generator should be of order 1 or 2 */
        nh = mat_mul(K[0][anz[0]-1],K[0][anz[0]-1]);
        Check_mat(nh);
        if (!nh->flags.Symmetric ||
            !nh->flags.Diagonal ||
            nh->array.SZ[0][0] != 1){
          finite_flag = FALSE;
        }
        free_mat(nh);
  
        /* the new generator should commute with all the other generators.
           (it suffices to check only those which have smaller index)
           (bare in mind that all of these have order at most 2) */
        for (i=0;i<anz[0]-1 && finite_flag;i++){
           nh = mat_kon(K[0][anz[0]-1],K[0][i],K[0][anz[0]-1]);
           if (mat_comp(nh,K[0][i]) != 0){
              finite_flag = FALSE;
           }
           free_mat(nh);
        }

        if (finite_flag == FALSE){
           return -1;
        }

      }
      else{
        free_mat(h);
      }

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

         if (mat_comp(h,erg[l]->representatives[tmp2])!=0){
            /* hinv = mat_inv(h);
            mat_muleq(hinv,erg[l]->representatives[tmp2]);
            free_mat(h);
            h = hinv; */
            hinv = mat_mul(erg[l]->rep_invs[tmp2],h);
            free_mat(h);
            h = hinv;
            if (l+1 < tmp){
               new_vec = mat_mul(h,erg[l+1]->orbit[0]);
               modp_mat(new_vec,2);
            }
         }
         else{
            /* only clean up memory */
            free_mat(h);
            return FALSE;
         }

         /* enter the next stage */
         flag = einordnen_2(erg,h,new_vec,++l,TRUE,K,anz);
      }
      else{
         if (l != (-1)){
            if (erg[l]->speicher<=erg[l]->length){
               extend_bahn(erg+l);
            }
            erg[l]->orbit[erg[l]->length]=new_vec;
            erg[l]->representatives[erg[l]->length]=h;
            /* inserted to reduce the number of mat_inv's */
            erg[l]->rep_invs[erg[l]->length] = mat_inv(h);
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

   if (flag == (-1)){
      return (-1);
   }

   if (flag){
      /* do a new orbit-calculation */
      /* apply all generators to the previous orbit elements */
      for(tmp=0;tmp<erg[l]->length;tmp++){
         for (tmp2=l;tmp2<erg[0]->representatives[0]->cols;tmp2++){
            for (i=erg[tmp2]->gen_no-1;i>=0;i--){
               h = erg[tmp2]->generators[i];
               nv = mat_mul(h,erg[l]->orbit[tmp]);
               modp_mat(nv,2);
               nh = mat_mul(h,erg[l]->representatives[tmp]);
               if (einordnen_2(erg,nh,nv,l,FALSE,K,anz) == -1){
                  return -1;
               }
            }
         }
      }
   }

   return flag;
} /* einordnen_2(..........) */

/***********************************************************************
@
@-----------------------------------------------------------------------
@
@ int strong_generators_2(matrix_TYP **base,bravais_TYP *U,
@                            matrix_TYP ***K,int *anz,MP_INT *mp)
@
@  This function decides whether or not the integral group generated
@  by the matrices in U->gen[i], 0<= i< U->gen_no, is finite.
@  If so, it will return TRUE and the order of the group via the multiple
@  precision integer mp. Otherwise the function will return FALSE immediately
@  if infiniteness can be proven.
@  In the case the given group is finite it will return a generating set
@  for the kernel of the Minkowski homomorphim (g -> g mod 2) via K,
@  and the number of matrices via anz[0].
@  In case the function proves the group to be infinite, it will return
@  at least one element in K, which proves the group to ge infinite.
@
@  matrix_TYP **base: a set of U->dim integral column vectors. there are
@                     two conditions to these vectors, which are not checked:
@                     1. the vectors entries 0,1 only,
@                     2. the vectors form a basis of Z^U->dim/(2Z)^U->dim.
@  bravais_TYP *U   : the group in question. It is in no way changed.
@  matrix_TYP ***K  : the function returns elements of the minkowski kernel
@                     via this pointer.
@  int *anz         : anz[0] is the number of returned elements in K
@  MP_INT *mp       : a multiple precesion integer. Call the function with
@                     strong_generators_2(....,&x) for an MP_INT x, which
@                     has been initialized via mpz_init(&x) before.
@
@ SIDEEFECT: none are known.
@-----------------------------------------------------------------------
@
************************************************************************/
int strong_generators_2(matrix_TYP **base,bravais_TYP *U,
                             matrix_TYP ***K,int *anz,MP_INT *mp)
{
   bahn **erg;

   int i,
       j,
       k,
       m,
       dim = U->dim,
       finite_flag,       /* indicates whether the kernel K stays finite */
       opt[6];            /* options for orbit_alg */

   bravais_TYP G;

   matrix_TYP **gen=U->gen,
              **con,
               *new_vec,
               *g,
               *h,
               *tmp;


   /* geting the memory for the result */
   erg = (bahn **) calloc(dim,sizeof(bahn *));
   for (i=0;i<dim;i++){
      erg[i] = (bahn *) calloc(1,sizeof(bahn));
      init_bahn(erg[i]);
   }

   /* initialise the minkowski kernel: to avoid trouble 1 is allways
      a generator */
   K[0] = (matrix_TYP **) malloc(MIN_SPEICHER*sizeof(matrix_TYP *));
   anz[0] = 0;
   finite_flag = TRUE;

   for (i=0;i<dim;i++){
      erg[i]->representatives[0] = init_mat(U->dim,U->dim,"1");
      erg[i]->rep_invs[0] = init_mat(U->dim,U->dim,"1");
      erg[i]->orbit[0] = copy_mat(base[i]);
      erg[i]->length = 1;
      erg[i]->hash->no = 0;
   }

   /* doing the orbit-calculation */
   for (j=0;j<erg[0]->length && finite_flag;j++){
     for (k=0;k<U->gen_no && finite_flag;k++){

        if (INFO_LEVEL & 4){
        /* printing all results for debugging */
        for (i=0;i<dim;i++){
           for (m=0;m<erg[i]->length;m++){
              printf("i=%d m=%d\n",i,m);
              if (erg[i]->representatives[m]->rows != dim ||
                  erg[i]->orbit[m]->rows != dim){
                printf("fehler in strong_generators_2\n");
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
        modp_mat(new_vec,2);

        //int old_anz = anz[0];
        if (einordnen_2(erg,h,new_vec,0,FALSE,K,anz) == -1){
            finite_flag = FALSE;
        }

        /* we got new generators for the kernel, so see
           whether we can prove the group to be infinite
        for (i=old_anz;i<anz[0] && finite_flag;i++){ */

           /* the new generator should be of order 1 or 2
           tmp = mat_mul(K[0][i],K[0][i]);
           Check_mat(tmp);
           if (!tmp->flags.Symmetric ||
               !tmp->flags.Diagonal ||
               tmp->array.SZ[0][0] != 1){
             finite_flag = FALSE;
           }
           free_mat(tmp); */

           /* the new generator should commute with all the other generators.
              (it suffices to check only those which have smaller index)
              (bare in mind that all of these have order at most 2)
           for (m=0;m<i && finite_flag;m++){
              tmp = mat_kon(K[0][i],K[0][m],K[0][i]);
              if (mat_comp(tmp,K[0][m]) != 0){
                 finite_flag = FALSE;
              }
              free_mat(tmp);
           }
        } */
      }
    }

   /* if the kernel K is trivial, we won't have generators at all,
      which causes trouble afterwards */
   if (anz[0] == 0){
      anz[0] = 1;
      K[0][0] = init_mat(U->dim,U->dim,"1");
   }

   /* this is the time to calculate the order of the group modulo
      the minkowski kernel */
   mpz_set_si(mp,1);
   for (i=0;i<dim;i++){
     mpz_mul_ui(mp,mp,(unsigned long) erg[i]->length);
   }

   /* bare in mind to free erg, which sounds crazy but is this
      way because the function was used in a different way before */
   for (i=0;i<dim;i++){
      free_bahn(erg[i]);
      free(erg[i]);
   }
   free(erg);

   /* we now decided whether nor not the group is finite,
      because if we still got the finite_flag, the kernel K
      is finitely generated, abelian, of exponent 2.
      bare in mind that we only got an generating system for
      K which generates K as normal subgroup of U */
   opt[0] = 4; opt[2] = 1; opt[1] = opt[3] = opt[4] = opt[5] = 0;
   for (i=0;i<dim;i++) opt[2] *= 2;
   if (finite_flag){
      /* make the normal closure */
      int old_anz = anz[0];
      for (i=0;i<old_anz && finite_flag;i++){
         con = orbit_alg(K[0][i],U,NULL,opt,&k);
         for (j=0;j<k;j++){
            for (m=0;m<anz[0] && mat_comp(K[0][m],con[j]) != 0;m++);
            if (m==anz[0]){
               if ((anz[0] % MIN_SPEICHER == 0) && anz[0] != 0)
                  K[0] = (matrix_TYP **) realloc(K[0],
                           (anz[0]+MIN_SPEICHER) * sizeof(matrix_TYP *));
               K[0][anz[0]] = con[j];
               Check_mat(K[0][anz[0]]);
               /* look whether it stays finite */
               if (finite_flag){
                  /* the new generator is of order 1 or 2 because it's
                     a conjugate of such an element, and the new generator
                     should commute with all the other generators. */
                  for (m=0;m<anz[0] && finite_flag;m++){
                     tmp = mat_kon(K[0][anz[0]],K[0][m],K[0][anz[0]]);
                     if (mat_comp(tmp,K[0][m]) != 0){
                        finite_flag = FALSE;
                     }
                     free_mat(tmp);
                  }
               }
               anz[0]++;
            }
            else{
               free_mat(con[j]);
            }

         }
         free(con);
      }
   }

   /* the finite_flag might be changed since we last checked it */   
   if (finite_flag){
      erg = NULL;
      G.gen = K[0];
      G.gen_no = anz[0];
      G.dim = dim;
      i = red_gen(&G,base,&erg,0);
      anz[0] = G.gen_no;
      K[0] = G.gen;
      mpz_mul_ui(mp,mp,(unsigned long) i);
      for (i=0;i<dim;i++){
         free_bahn(erg[i]);
         free(erg[i]);
      }
      free(erg);
   }

   return finite_flag;
} /* strong_generators_2 */
