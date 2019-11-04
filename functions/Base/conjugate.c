#include"typedef.h"
#include"longtools.h"
#include"matrix.h"
#include"base.h"

/**************************************************************************
@
@ -------------------------------------------------------------------------
@ FILE: conjugate.c
@ -------------------------------------------------------------------------
@
***************************************************************************/

/***********************************************************************
@
@ -------------------------------------------------------------------------
@
@ matrix_TYP *conjugated(bravais_TYP *G,bravais_TYP *H,
@                        matrix_TYP **N,int Nanz,bahn **strong)
@ Suppose the matrices in G->gen, H->gen generate finite matrix groups, and N
@ generates a matrix group with the condition that the orbit of N
@ on the conjugated of G under N is finite. Then the algorithm decides
@ whether G is conjugated to a supergroup of H, and if so gives an element
@ conjugating G over H, i.e H<= X * G * X^-1 for the result of this function.
@ Declaration of variables:
@      G:       matrices generating G,
@      H:
@      N:
@      strong:   a set of base/strong generators for G returned by
@                strong_generators
@
@ Sideefects: None of the variables given nor global variables should
@             be affected.
@
@ -------------------------------------------------------------------------
@
***********************************************************************/ 
matrix_TYP *conjugated(bravais_TYP *G,bravais_TYP *H,
                       matrix_TYP **N,int Nanz,bahn **strong)
{

matrix_TYP **orbit,     /* represents the orbit of G under N by the
                           conjugating element */
            *test_subgroup,
            *tmp,
            *tmp2,
            *tmp3,
            *tmp4,
            *erg = NULL;

int found=1,  /* number of conjugated subgroups found so far */
    dealt=0,  /* number of those dealt with */
    i,
    j,
    k,
    flag;

   if (INFO_LEVEL & 4){
      fprintf(stderr,"entering conjugated\n");
   }


   /* let's see whether H is not already a Subgroup of G */
   flag = TRUE;
   k = 0;
   while ((k<H->gen_no) && flag){
      if (!is_element(H->gen[k],G,strong,NULL)){
         flag = FALSE;
      }
      k++;
   }
   /* if it was, leave instantly */
   if (flag){
      erg = init_mat(G->gen[0]->cols,G->gen[0]->cols,"1");
      return(erg);
   }

   orbit = (matrix_TYP **) malloc(1*sizeof(matrix_TYP *));
   orbit[0] = init_mat(G->gen[0]->cols,G->gen[0]->cols,"1");

   /* start the obit calculation */
   while ((dealt<found) && (erg == NULL)){

     /* output for debugging */
     if (INFO_LEVEL & 4){
        fprintf(stderr,"got %i conjugated subgroups so far\n",found);
        fprintf(stderr,"dealt with %i of these\n",dealt);
     }
     /* loop over the generators of N */
     for (i=0;(i<Nanz) & (erg == NULL);i++){

        /* calculating the new element */
        test_subgroup = mat_mul(orbit[dealt],N[i]);

        /* see if this gives really a new element */
        j=0;
        tmp = mat_inv(test_subgroup);
        while ((test_subgroup != NULL) && (j< found)){
           tmp2 = mat_mul(orbit[j],tmp);
           tmp3 = mat_inv(tmp2);
           k=0;
           flag = TRUE;
           while (flag && (k<G->gen_no)){
              tmp4 = mat_mul(tmp2,G->gen[k]);
              tmp4 = mat_muleq(tmp4,tmp3);
              if (!is_element(tmp4,G,strong,NULL)){
                 flag = FALSE;
              }
              free_mat(tmp4);
              k++;
           }

           if (flag){
              free_mat(test_subgroup);
              test_subgroup = NULL;
           }
           free_mat(tmp2);
           free_mat(tmp3);
           j++;
        }

        /* if it was a new element, we should memorize it and check whether
           is over H */
        if ( test_subgroup!= NULL){
           /* so let's see whether it is really over H */
           flag = TRUE;
           k=0;
           while(flag && (k<H->gen_no)){
              tmp2 = mat_mul(test_subgroup,H->gen[k]);
              mat_muleq(tmp2,tmp);
              if (!is_element(tmp2,G,strong,NULL)){
                 flag = FALSE;
              }
              free_mat(tmp2);
              k++;
           }
           /* get more memory */
           orbit = (matrix_TYP **) realloc(orbit,(found+1)*sizeof(matrix_TYP*));
           orbit[found] = test_subgroup;
           found++;
           if (flag){ erg = long_mat_inv(test_subgroup); }
        }

        free_mat(tmp);

     } /* for (i= .... */
     dealt++;
   }


   /* cleaning up memory */
   for (i=0;i<found;i++){
      free_mat(orbit[i]);
   }
   free(orbit);

   return erg;
}

