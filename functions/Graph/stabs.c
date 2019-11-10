/* last change: 28.09.00 by Oliver Heidbuechel */


#include<typedef.h>
#include<getput.h>
#include<matrix.h>
#include<longtools.h>
#include<tools.h>
#include"zass.h"
#include <bravais.h>
#include <sort.h>
#include <presentation.h>
#include <base.h>
#include <graph.h>


boolean GRAPH_DEBUG;

/* -------------------------------------------------------------------------- */
/* returns 1 iff mat = Id                                                     */
/* -------------------------------------------------------------------------- */
static int is_id(matrix_TYP *mat)
{
   int i, j;

   for (i = 0; i < mat->rows; i++){
      for (j = 0; j < mat->cols; j++){
         if (i == j && mat->array.SZ[i][i] != 1)
            return(0);
         if (i != j && mat->array.SZ[i][j] != 0)
            return(0);
      }
   }
   return(1);
}



/* -------------------------------------------------------------------- */
/* Setze Matrizen (gegeben in N) in Worte ein.                          */
/* Streiche Identitaet und doppelte.                                    */
/* -------------------------------------------------------------------- */
/* words: Worte                                                         */
/* no_of_words: Anzahl der Worte                                        */
/* N: Matrizen, die eingesetzt werden sollen                            */
/*    (muss zu words "passen")                                          */
/* N_inv: Inversen von N (werden bei Bedarf berechnet)                  */
/* anz: Anzahl von N                                                    */
/* flag: wenn = 0, wird einfach N zurueckgegeben                        */
/* stab_gen_no: hier wird die Anzahl der Elemente zurueckgegeben        */
/* -------------------------------------------------------------------- */
matrix_TYP **stab_coz(int **words,
                      int no_of_words,
                      matrix_TYP **N,
                      matrix_TYP **N_inv,
                      int anz,
                      int flag,
                      int *stab_gen_no)
{
   int i, counter = 0;

   matrix_TYP **stab;



   if (flag == 0){
      /* the stabilizer of the cocycle is the full normalizer */
      stab = (matrix_TYP **)calloc(anz, sizeof(matrix_TYP *));
      stab_gen_no[0] = anz;
      for (i = 0; i < anz; i++){
         stab[i] = copy_mat(N[i]);
      }
   }
   else{
      stab = (matrix_TYP **)calloc(no_of_words, sizeof(matrix_TYP *));
      for (i = 0; i < no_of_words; i++){
         normalize_word(words[i]);
         stab[counter] = mapped_word(words[i], N, N_inv);

         if (!is_id(stab[counter])){
            if (mat_search(stab[counter], stab, counter, mat_comp) != -1){
               free_mat(stab[counter]);
            }
            else{
               mat_quicksort(stab, 0, counter, mat_comp);
               counter++;
            }
         }
         else{
            free_mat(stab[counter]);
         }
      }
      stab_gen_no[0] = counter;
   }

   return(stab);
}



/* -------------------------------------------------------------------------- */
/* Did we have this lattice bevore? If not, insert the lattice in the tree!   */
/* -------------------------------------------------------------------------- */
static int search_in_tree(struct tree *Baum,
                          matrix_TYP **list,
                          matrix_TYP *lattice,
                          int number)
{
   int flag, no;

   no = Baum->no;
   flag = cmp_mat(list[no], lattice);

   if (flag == 0)
      return(no);
   if (flag < 0){
      if (Baum->left == NULL){
         Baum->left = (struct tree *)calloc(1, sizeof(struct tree));
         Baum->left->no = number;
         return(-1);
      }
      else{
         return(search_in_tree(Baum->left, list, lattice, number));
      }
   }
   else{
      if (Baum->right == NULL){
         Baum->right = (struct tree *)calloc(1, sizeof(struct tree));
         Baum->right->no = number;
         return(-1);
      }
      else{
         return(search_in_tree(Baum->right, list, lattice, number));
      }

   }
}



/* -------------------------------------------------------------------------- */
/* mark the lattices in the orbit under the stablizer of a cocycle            */
/* -------------------------------------------------------------------------- */
static void mark_lattice(boolean *lattice_orbit,
                         matrix_TYP **LIST,
                         int anz,
                         matrix_TYP *mat)
{
   int i;

   for (i = 0; i < anz; i++){
      if (cmp_mat(mat, LIST[i]) == 0){
         lattice_orbit[i] = TRUE;
         return;
      }
   }
   fprintf(stderr, "ERROR in mark_lattice!\n");
   exit(7);
}



/* -------------------------------------------------------------------------- */
/* Calculate the stabilizer of a lattice and return the words;                */
/* save the number of words in no                                             */
/* -------------------------------------------------------------------------- */
int **stab_lattice(matrix_TYP *lattice,
                   matrix_TYP **N,
                   int N_gen_no,
                   int *no,
                   matrix_TYP **LIST,
                   int anz,
                   boolean *lattice_orbit)
{
   int i, j, k, flag,
       counter = 1, number,
       **words, **WORDS;

   matrix_TYP **list;

   struct tree *Baum;


   list = (matrix_TYP **)calloc(MIN_SPEICHER, sizeof(matrix_TYP *));
   words = (int **)calloc(MIN_SPEICHER, sizeof(int *));
   WORDS = (int **)calloc(MIN_SPEICHER, sizeof(int *));
   list[0] = lattice;
   Baum = (struct tree *)calloc(1, sizeof(struct tree));
   Baum->no = i = no[0] = 0;

   while (i < counter){
      for (j = 0; j < N_gen_no; j++){
         list[counter] = mat_mul(N[j], list[i]);
         long_col_hnf(list[counter]);
         flag = search_in_tree(Baum, list, list[counter], counter);
         if (flag == -1){
            /* new lattice */
            words[counter] = (int *)calloc(MIN_SPEICHER, sizeof(int));
            if (i != 0){
               memcpy(words[counter], words[i], (words[i][0] + 1) * sizeof(int));
               words[counter][words[i][0] + 1] = -(j+1);
               words[counter][0]++;
            }
            else{
               words[counter][0] = 1;
               words[counter][1] = -(j+1);
            }
            if (lattice_orbit != NULL && LIST != NULL && anz != 0)
               mark_lattice(lattice_orbit, LIST, anz, list[counter]);
            counter++;
         }
         else{
            /* we got this lattice already */
            free_mat(list[counter]);
            WORDS[no[0]] = (int *)calloc(MIN_SPEICHER, sizeof(int));
            if (i != 0){
               memcpy(WORDS[no[0]], words[i], (words[i][0] + 1) * sizeof(int));
               WORDS[no[0]][words[i][0] + 1] = -(j+1);
               WORDS[no[0]][0]++;
            }
            else{
               WORDS[no[0]][0] = 1;
               WORDS[no[0]][1] = -(j+1);
            }
            if (flag != 0){
               number = WORDS[no[0]][0] + words[flag][0];
               for (k = 0; k < words[flag][0]; k++){
                  WORDS[no[0]][number - k] = -words[flag][k + 1];
               }
               WORDS[no[0]][0] = number;
            }
            no[0]++;
         }
      }
      i++;
   }

   /* clean */
   for (i = 1; i < counter; i++){
      free(words[i]);
      free_mat(list[i]);
   }
   free(words);
   free(list);
   free_tree(Baum);

   return(WORDS);
}



/* -------------------------------------------------------------------------- */
/* returns 1 if  word[n] == word[i] for a 0 <= i < n                          */
/* returns 0 otherwise                                                        */
/* -------------------------------------------------------------------------- */
int word_already_there(int **WORDS,
                       int n)
{
   int i, j, flagge = 0;

   for (i = 0; i < n && flagge == 0; i++){
      flagge = 1;
      for (j = 0; j <= WORDS[i][0]; j++){
         if (WORDS[i][j] != WORDS[n][j]){
            flagge = 0;
            break;
         }
      }
   }
   return(flagge);
}



/* -------------------------------------------------------------------------- */
/* Calculate the stabilizer of lattice under N and return the words;          */
/* save the number of words in no                                             */
/* -------------------------------------------------------------------------- */
matrix_TYP **calculate_S1(matrix_TYP *lattice,
                          matrix_TYP **N,
                          int anz,
                          int *no,
                          matrix_TYP **dataN,
                          matrix_TYP **dataNinv,
                          matrix_TYP *dataX)
{
   int i, j, k, flag, i__, addmem,
       counter = 1, number,
       **words, **WORDS;

   matrix_TYP **list, **S1,
              **Ninv, *tmp;

   struct tree *Baum;

   if (GRAPH_DEBUG){
      Ninv = (matrix_TYP **)calloc(anz, sizeof(matrix_TYP *));
   } else
      Ninv = 0;

   list = (matrix_TYP **)calloc(1024, sizeof(matrix_TYP *));
   words = (int **)calloc(1024, sizeof(int *));
   WORDS = (int **)calloc(1024, sizeof(int *));
   list[0] = lattice;
   Baum = (struct tree *)calloc(1, sizeof(struct tree));
   Baum->no = i = no[0] = 0;

   while (i < counter){
      for (j = 0; j < anz; j++){
         if (counter % 1024 == 0){
	    list = (matrix_TYP **)realloc(list, (counter + 1024) * sizeof(matrix_TYP *));
	 }
         list[counter] = mat_mul(N[j], list[i]);
         long_col_hnf(list[counter]);
         flag = search_in_tree(Baum, list, list[counter], counter);
         if (flag == -1){
            /* new lattice */
	    if (counter % 1024 == 0){
	       words = (int **)realloc(words, (counter + 1024) * sizeof(int *));
	    }
            /* words[counter] = (int *)calloc(1024, sizeof(int)); */
            if (i != 0){
	       words[counter] = (int *)calloc(words[i][0] + 2, sizeof(int));
               /*memcpy(words[counter], words[i], (words[i][0] + 1) * sizeof(int));*/
	       for (i__ = 0; i__ < words[i][0] + 1; i__++){
	          words[counter][i__] = words[i][i__];
	       }
               words[counter][words[i][0] + 1] = -(j+1);
               words[counter][0]++;
            }
            else{
	       words[counter] = (int *)calloc(2, sizeof(int));
               words[counter][0] = 1;
               words[counter][1] = -(j+1);
            }
            counter++;
         }
         else{
            /* we got this lattice already */
            free_mat(list[counter]);
	    if (no[0] % 1024){
	       WORDS = (int **)realloc(WORDS, (no[0] + 1024) * sizeof(int *));
	    }
            /*WORDS[no[0]] = (int *)calloc(1024, sizeof(int));*/
	    if (flag != 0){
	       addmem = words[flag][0];
	    }
	    else{
	       addmem = 0;
	    }
            if (i != 0){
               WORDS[no[0]] = (int *)calloc(words[i][0] + 2 + addmem, sizeof(int));
               /*memcpy(WORDS[no[0]], words[i], (words[i][0] + 1) * sizeof(int));*/
	       for (i__ = 0; i__ < words[i][0] + 1; i__++){
	          WORDS[no[0]][i__] = words[i][i__];
	       }
               WORDS[no[0]][words[i][0] + 1] = -(j+1);
               WORDS[no[0]][0]++;
            }
            else{
	       WORDS[no[0]] = (int *)calloc(2 + addmem, sizeof(int));
               WORDS[no[0]][0] = 1;
               WORDS[no[0]][1] = -(j+1);
            }
            if (flag != 0){
               number = WORDS[no[0]][0] + words[flag][0];
               for (k = 0; k < words[flag][0]; k++){
                  WORDS[no[0]][number - k] = -words[flag][k + 1];
               }
               WORDS[no[0]][0] = number;
            }
            no[0]++;
         }
      }
      i++;
   }

   if (GRAPH_DEBUG){
      /* is this really an element of the stabilizer? */
      for (i = 0; i < no[0]; i++){
         tmp = mapped_word(WORDS[i], N, Ninv);
         tmp = mat_muleq(tmp, lattice);
         long_col_hnf(tmp);
         if (cmp_mat(tmp, lattice) != 0){
            fprintf(stderr, "MIST!\n");
            exit(9);
         }
         free_mat(tmp);
      }
      for (i = 0; i < anz; i++){
         if (Ninv[i] != NULL) free_mat(Ninv[i]);
      }
      free(Ninv);
   }

   /* clean */
   for (i = 1; i < counter; i++){
      free(words[i]);
      free_mat(list[i]);
   }
   free(words);
   free(list);
   free_tree(Baum);

   /* now calculate the matrices */
   counter = 0;
   S1 = (matrix_TYP **)calloc(no[0], sizeof(matrix_TYP *));
   for (i = 0; i < no[0]; i++){
      /* simplify the words */
      normalize_word(WORDS[i]);
      if (WORDS[i][0] == 1 && WORDS[i][1] < 0)
         WORDS[i][1] *= (-1);
      if (word_already_there(WORDS, i) == 0){
         S1[counter] = graph_mapped_word(WORDS[i], dataN, dataNinv, dataX);
         if (mat_search(S1[counter], S1, counter, mat_comp) != -1){
            free_mat(S1[counter]);
         }
         else{
            mat_quicksort(S1, 0, counter, mat_comp);
            counter++;
         }
      }
   }
   for (i = 0; i < no[0]; i++)
      free(WORDS[i]);
   free(WORDS);

   no[0] = counter;
   return(S1);
}





