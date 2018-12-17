/* last change: 01.03.01 by Oliver Heidbuechel */


#include <typedef.h>
#include <matrix.h>
#include <bravais.h>
#include <base.h>
#include <graph.h>
#include <presentation.h>
#include <longtools.h>
#include <datei.h>


extern int INFO_LEVEL;


/* -------------------------------------------------------------------- */
/* transform the i-th row of a matrix to a list of integers             */
/* -------------------------------------------------------------------- */
static int *row_to_list(matrix_TYP *mat,
                        int i)
{
   int *list, n;

   list = (int *)calloc(mat->cols, sizeof(int));
   for (n = 0; n < mat->cols; n++)
      list[n] = mat->array.SZ[i][n];

   return(list);
}



/* -------------------------------------------------------------------- */
/* in mats we have words for the generators of one represenative of     */
/* each conjugacy class of maximal subgroups of G                       */
/* in the last row there has to be the number of elements in the        */
/* conjugacy class                                                      */
/* fill a t_sub_TYP - structure with the generators for each            */
/* representative of a conjugacy class                                  */
/* -------------------------------------------------------------------- */
static t_sub_TYP *t_sub_info(matrix_TYP **mats,
                             int anz,
                             bravais_TYP *G,
                             matrix_TYP **G_geninv)
{
   t_sub_TYP *sub;

   int i, j, no;

   matrix_TYP **base;


   /* prepare */
   sub = (t_sub_TYP *)calloc(1, sizeof(t_sub_TYP));
   sub->con_no = anz;
   sub->groups = (bravais_TYP **)calloc(anz, sizeof(bravais_TYP *));
   sub->worte = (int ***)calloc(anz, sizeof(int **));
   sub->elem_no = (int *)calloc(anz, sizeof(int));
   sub->strong = (bahn ***)calloc(anz, sizeof(bahn **));
   sub->orders = (int *)calloc(anz, sizeof(int));
   sub->total = 0;

   for (i = 0; i < anz; i++){
      no = mats[i]->rows - 1;
      sub->groups[i] = init_bravais(G->dim);
      if (no == 0)
         sub->groups[i]->gen = NULL;
      else
         sub->groups[i]->gen = (matrix_TYP **)calloc(no, sizeof(matrix_TYP *));
      sub->groups[i]->gen_no = no;
      sub->worte[i] = (int **)calloc(no, sizeof(int *));
      for (j = 0; j < no; j++){
         sub->worte[i][j] = row_to_list(mats[i], j);
         sub->groups[i]->gen[j] = mapped_word(sub->worte[i][j], G->gen, G_geninv);
      }
      /* in the last row there is the number of elements in the conjugacy class
         and the order of the elements */
      sub->elem_no[i] = mats[i]->array.SZ[no][0];
      sub->orders[i] = sub->groups[i]->order = mats[i]->array.SZ[no][1];
      sub->total += sub->elem_no[i];

      /* calculate strong generating set */
      base = get_base(G);
/*
???????????????????????????????????????????????????????????????????????
*/
      sub->strong[i] = strong_generators(base, sub->groups[i], TRUE);
      for (j = 0; j < G->dim; j++){
         free_mat(base[j]);
      }
      free(base);
   }

   return(sub);
}



/* -------------------------------------------------------------------- */
/* free a t_sub_TYP - structure                                         */
/* -------------------------------------------------------------------- */
static void free_t_sub_TYP(t_sub_TYP *sub)
{
   int i, j, dim = sub->groups[0]->dim;


   for (i = 0; i < sub->con_no; i++){
      for (j = 0; j < sub->groups[i]->gen_no; j++){
         free(sub->worte[i][j]);
      }
      free(sub->worte[i]);
      free_bravais(sub->groups[i]);
      for (j = 0; j < dim; j++){
         free_bahn(sub->strong[i][j]);
         free(sub->strong[i][j]);
      }
      free(sub->strong[i]);
   }
   free(sub->strong);
   free(sub->worte);
   free(sub->groups);
   free(sub->elem_no);
   free(sub->orders);
   free(sub);
}



/* --------------------------------------------------------------------
Suppose the matrices in G->gen generate a finite matrix group, and N
generates a matrix group with the condition that the orbit of N
on the conjugated of G under N is finite.

The groups in the orbit will be given by the element that conjugates them!

Declaration of variables:
   G:       matrices generating G,
   N:       see above
   Nanz:    number of elements in N
   strong:  a set of base/strong generators for G returned by strong_generators

Sideefects: None of the variables given nor global variables should be affected.

The code is out of Base/conjugate.c!
-------------------------------------------------------------------- */
static matrix_TYP **conjugated_groups(bravais_TYP *G,
                                      matrix_TYP **N,
                                      int Nanz,
                                      bahn **strong,
                                      int *anz)
{

matrix_TYP **orbit,     /* represents the orbit of G under N by the
                           conjugating element */
            *test_subgroup,
            *tmp,
            *tmp2,
            *tmp3,
            *tmp4;

int found=1,  /* number of conjugated subgroups found so far */
    dealt=0,  /* number of those dealt with */
    i,
    j,
    k,
    flag;

   if (INFO_LEVEL & 4){
      fprintf(stderr,"entering conjugated\n");
   }


   orbit = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
   orbit[0] = init_mat(G->gen[0]->cols, G->gen[0]->cols, "1");

   /* start the orbit calculation */
   while (dealt < found){

     /* output for debugging */
     if (INFO_LEVEL & 4){
        fprintf(stderr,"got %i conjugated subgroups so far\n",found);
        fprintf(stderr,"dealt with %i of these\n",dealt);
     }
     /* loop over the generators of N */
     for (i=0; i < Nanz;i++){

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

        /* get more memory */
        if (test_subgroup != NULL){
           orbit = (matrix_TYP **) realloc(orbit, (found + 1) * sizeof(matrix_TYP*));
           orbit[found] = test_subgroup;
           found++;
        }
        free_mat(tmp);

     } /* for (i= .... */
     dealt++;
   }

   anz[0] = found;

   return(orbit);
}



/* -------------------------------------------------------------------- */
/* Let G and H be two groups with the same finite order.                */
/* strong has to be the result of a call of strong_generators() for G   */
/* returns TRUE iff H == G                                              */
/* -------------------------------------------------------------------- */
static boolean equal_groups(bravais_TYP *H,
                            bravais_TYP *G,
                            bahn **strong)
{
   int i;

   boolean FLAG = TRUE;


   for (i = 0; FLAG && i < H->gen_no; i++){
      if (!is_element(H->gen[i], G, strong, NULL)) FLAG = FALSE;
   }

   return(FLAG);
}



/* -------------------------------------------------------------------- */
/* If the elements in conjugacy class i and j are in one orbit, the     */
/* conjugacy classes have to have the same number of elements and the   */
/* representatives have to have the same order.                         */
/* -------------------------------------------------------------------- */
static int *construct_help_list(t_sub_TYP *sub)
{
   int *help, i, j, counter = 1;

   boolean flag;


   help = (int *)calloc(sub->con_no, sizeof(int));
   for (i = 0; i < sub->con_no; i++){
      if (help[i] == 0){
         help[i] = -counter;
         flag = FALSE;
         for (j = i + 1; j < sub->con_no; j++){
            if (sub->orders[j] == sub->orders[i] &&
                sub->elem_no[j] == sub->elem_no[i]){
               help[j] = -counter;
               flag = TRUE;
            }
         }
         if (flag == FALSE){
            /* there is only one conjugacy class with this combination
               of order and length */
            help[i] = 0;
         }
         counter++;
      }
   }

   return(help);
}



/* -------------------------------------------------------------------- */
/* mats: matrices with information about the maximal subgroups of G out */
/*       of a GAP-programm                                              */
/* no: no of mats                                                       */
/* pres: presentation for G                                             */
/* aff_no: save the number of affine classes here                       */
/* aff_cons: save the no of conjugacy classes of each affine class here */
/* aff_cons_no: save the no of elements in each conjugacy class of each */
/*              affine class here                                       */
/* R: save representatives for the affine classes here                  */
/* flag: 0: only representatives are calculated                         */
/*       1: all subgroups are calculated                                */

/* G->normal and G->gen have to generate the normalizer of G!!!!!       */

/* -------------------------------------------------------------------- */
bravais_TYP ****t_subgroups(bravais_TYP *G,
                            matrix_TYP **mats,
                            int no,
                            matrix_TYP *pres,
                            int *aff_no,
                            int **aff_cons,
                            int ***aff_cons_no,
                            bravais_TYP ***R,
                            int flag)
{
   t_sub_TYP *sub;

   bravais_TYP ****S, **group;

   int i, j, k, l, counter, together,
       cen_no, orbit_no, coho_size,
       *stab_genno, ***WORDS, *NUMBER_OF_WORDS, *trashlist,
       *help_list;

   matrix_TYP **G_geninv, **coz, **X, **trash, **G_norminv, ***stab,
              **conj, **Rinv, *inv, *cocycle;

   MP_INT *names;


   /* deal with trivial case */
   if (no == 0){
      /* there is no subgroup */
      R[0] = (bravais_TYP **)calloc(1, sizeof(bravais_TYP *));
      R[0][0] = init_bravais(G->dim + 1);
      R[0][0]->gen = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
      R[0][0]->gen_no = 1;
      R[0][0]->gen[0] = init_mat(G->dim + 1, G->dim + 1, "1");
      aff_cons[0] = NULL;
      aff_cons_no[0] = NULL;
      aff_no[0] = 0;
      return(NULL);
   }

   /* prepare */
   cen_no = G->cen_no; G->cen_no = 0;
   G_geninv = (matrix_TYP **)calloc(G->gen_no, sizeof(matrix_TYP *));
   sub = t_sub_info(mats, no, G, G_geninv);
   WORDS = (int ***)calloc(MIN_SPEICHER, sizeof(int **));
   NUMBER_OF_WORDS = (int *)calloc(MIN_SPEICHER, sizeof(int));
   coz = all_cocycles(pres, G, aff_no, G_geninv, &X, &names, WORDS, NUMBER_OF_WORDS,
                      &trash, &coho_size, &trashlist, 1);
   if (trash != NULL){
      for (i = 0; i < G->normal_no; i++)
         if (trash[i] != NULL) free_mat(trash[i]);
      free(trash);
   }
   for (i = 0; i < G->gen_no; i++)
      if (G_geninv[i] != NULL) free_mat(G_geninv[i]);
   free(G_geninv);
   aff_cons[0] = (int *)calloc(aff_no[0], sizeof(int));
   aff_cons_no[0] = (int **)calloc(aff_no[0], sizeof(int *));
   R[0] = (bravais_TYP **)calloc(aff_no[0], sizeof(bravais_TYP *));
   for (i = 0; i < aff_no[0]; i++){
      aff_cons_no[0][i] = (int *)calloc(sub->con_no, sizeof(int));
      R[0][i] = extract_r(G, coz[i]);
   }
   S = (bravais_TYP ****)calloc(aff_no[0], sizeof(bravais_TYP ***));

   /* is it another trivial case */
   if (no == 1 && mats[0]->rows == 1){
      /* only the trivial group is a subgroup */
      for (i = 0; i < aff_no[0]; i++){
         S[i] = (bravais_TYP ***)calloc(1, sizeof(bravais_TYP **));
         S[i][0] = (bravais_TYP **)calloc(1, sizeof(bravais_TYP *));
         S[i][0][0] = init_bravais(G->dim + 1);
         S[i][0][0]->gen = (matrix_TYP **)calloc(1, sizeof(matrix_TYP *));
         S[i][0][0]->gen_no = 1;
         S[i][0][0]->gen[0] = init_mat(G->dim +  1, G->dim + 1, "1");
         aff_cons[0][i] = 1;
         aff_cons_no[0][i][0] = 1;
      }
   }
   else{
      /* get the point group of the affine normalizer for each affine class representative */
      G_norminv = (matrix_TYP **)calloc(G->normal_no, sizeof(matrix_TYP *));
      stab = (matrix_TYP ***)calloc(aff_no[0], sizeof(matrix_TYP **));
      stab_genno = (int *)calloc(aff_no[0], sizeof(int));
      for (i = 0; i < aff_no[0]; i++){
         stab[i] = stab_coz(WORDS[i], NUMBER_OF_WORDS[i], G->normal, G_norminv,
                            G->normal_no, i, &stab_genno[i]);
         stab[i] = realloc(stab[i], (stab_genno[i] + G->gen_no) * sizeof(matrix_TYP *));
         for (j = 0; j < G->gen_no; j++)
            stab[i][j + stab_genno[i]] = G->gen[j];
      }

      /* which subgroups are in the same orbit */
      for (i = 0; i < aff_no[0]; i++){
         Rinv = (matrix_TYP **)calloc(R[0][i]->gen_no, sizeof(matrix_TYP *));
         counter = 0;
         help_list = construct_help_list(sub);
         S[i] = (bravais_TYP ***)calloc(sub->con_no, sizeof(bravais_TYP **));
         for (j = 0; j < sub->con_no; j++){
            if (help_list[j] <= 0){

               if (help_list[j] < 0){
                  /* possibly we have to join this conjugacy class with another one */
                  help_list[j] *= -1;

                  /* calculate orbit */
                  conj = conjugated_groups(sub->groups[j], stab[i], stab_genno[i] + G->gen_no,
                                           sub->strong[j], &orbit_no);
                  group = (bravais_TYP **)calloc(orbit_no, sizeof(bravais_TYP *));

                  if (orbit_no % sub->elem_no[j] != 0){ /* paranoia test */
                     fprintf(stderr, "ERROR 1 in t_subgroups\n");
                     exit(6);
                  }

                  together = orbit_no / sub->elem_no[j];
                  /* we have to join "together" conjugacy classes - which? */
                  if (together > 1){
                     for (k = 1; together > 1 && k < orbit_no; k++){
                        inv = long_mat_inv(conj[k]);
                        group[k] = konj_bravais(sub->groups[j], inv);
                        free_mat(inv);
                        for (l = j + 1; together > 1 && l < sub->con_no; l++){
                           if (equal_groups(group[k], sub->groups[l], sub->strong[l])){
                              together--;
                              help_list[l] *= -1;
                           }
                        }
                     }
                  }
                  if (together != 1){	/* paranoia test */
                     fprintf(stderr, "ERROR 2 in t_subgroups\n");
                     exit(5);
                  }
               }
               else{
                  /* only the elements in this conjugacy class are in an orbit under
                     the normalizer */
                  if (sub->elem_no[j] == 1){
                     /* only one element in the conjugacy class */
                     conj = NULL;
                     orbit_no = 1;
                     group = NULL;
                  }
                  else{
                     conj = conjugated_groups(sub->groups[i], G->gen, G->gen_no,
                                              sub->strong[i], &orbit_no);
                     group = (bravais_TYP **)calloc(orbit_no, sizeof(bravais_TYP *));
                  }
               }
               if (flag == 0)
                  S[i][counter] = (bravais_TYP **)calloc(1, sizeof(bravais_TYP *));
               else
                  S[i][counter] = (bravais_TYP **)calloc(orbit_no, sizeof(bravais_TYP *));

               /* calculate the first element in the conjugacy class */
               S[i][counter][0] = init_bravais(R[0][i]->dim);
               S[i][counter][0]->gen_no = sub->groups[j]->gen_no;
               S[i][counter][0]->gen = (matrix_TYP **)calloc(sub->groups[j]->gen_no,
                                                             sizeof(matrix_TYP *));
               for (k = 0; k < sub->groups[j]->gen_no; k++){
                  S[i][counter][0]->gen[k] = mapped_word(sub->worte[j][k], R[0][i]->gen, Rinv);
               }

               if (flag != 0){
                  /* calculate the other groups in this conjugacy class */
                  for (k = 1; k < orbit_no; k++){
                     if (group[k] == NULL){
                        inv = long_mat_inv(conj[k]);
                        group[k] = konj_bravais(sub->groups[j], inv);
                        free_mat(inv);
                     }
                     cocycle = sg(R[0][i], group[k]);
                     S[i][counter][k] = extract_r(group[k], cocycle);
                     free_mat(cocycle);
                  }
               }
               aff_cons_no[0][i][counter] = orbit_no;
               counter++;

               for (k = 0; group != NULL && k < orbit_no; k++)
                  if (group[k] != NULL) free_bravais(group[k]);
               if (group != NULL) free(group);
               for (k = 0; conj != NULL && k < orbit_no; k++)
                  free_mat(conj[k]);
               if (conj != NULL) free(conj);
            }
         }
         free(help_list);
         aff_cons[0][i] = counter;
         for (j = 0; j < R[0][i]->gen_no; j++)
            if (Rinv[j] != NULL) free_mat(Rinv[j]);
         free(Rinv);
      }

      /* clean */
      for (i = 0; i < aff_no[0]; i++){
         for (j = 0; j < stab_genno[i]; j++)
            free_mat(stab[i][j]);
         free(stab[i]);
      }
      free(stab);
      free(stab_genno);
      for (i = 0; i < G->normal_no; i++)
         if (G_norminv[i] != NULL) free_mat(G_norminv[i]);
      free(G_norminv);
   }

   /* clean */
   free_t_sub_TYP(sub);
   for (i = 0; i < 3; i++)
      free_mat(X[i]);
   free(X);
   for (i = 0; i < aff_no[0]; i++){
      free_mat(coz[i]);
      if (WORDS[i] != NULL){
         for (j = 0; j < NUMBER_OF_WORDS[i]; j++)
            if (WORDS[i][j] != NULL)
               free(WORDS[i][j]);
         free(WORDS[i]);
      }
      mpz_clear(names + i);
   }
   free(coz);
   if (WORDS != NULL) free(WORDS);
   free(NUMBER_OF_WORDS);
   free(names);

   G->cen_no = cen_no;
   return(S);
}









