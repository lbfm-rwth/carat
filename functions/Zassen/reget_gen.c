#include "typedef.h"
#include "matrix.h"
#include "gmp.h"
#include "zass.h"
#include "getput.h"

extern int INFO_LEVEL;

static int and(int *a,int n)
{
  int i;

  for (i=0;i<n;i++) if (a[i] == FALSE) return FALSE;

  return TRUE;
}



static int is_equal_red(matrix_TYP *x,matrix_TYP *y)
{
  int i,
      j,
      rows,
      cols;

  if (x->cols < y->cols)
     cols = x->cols;
  else
     cols = y->cols;

  if (x->rows < y->rows)
     rows = x->rows;
  else
     rows = y->rows;

  for (i=0;i<cols;i++)
     for (j=0;j<rows;j++)
        if (x->array.SZ[i][j] != y->array.SZ[i][j]) return FALSE;

  return TRUE;
}


static int red_pos(matrix_TYP *x,matrix_TYP **A,int n)
{
  int i;

  for (i=0;i<n;i++)
    if (is_equal_red(x,A[i])) return i;

  return -1;

}

/**********************************************************************
@
@----------------------------------------------------------------------
@
@ matrix_TYP *reget_gen(matrix_TYP **map,int number,bravais_TYP *G,
@                              int **words,int word_flag)
@
@----------------------------------------------------------------------
@
***********************************************************************/
matrix_TYP *reget_gen(matrix_TYP **map,int number,bravais_TYP *G,
                             int **words,int word_flag)
{
   int *found,              /* gives a flag for each generator of G, TRUE iff
                               we already found this element */
        length = G->gen_no, /* the number of elements found so far */
        speicher = MIN_SPEICHER,
        k,
        i,
        j;

   int **ele_words;        /* for each matrix in ele we store a word in the
                              generators here which gives this element:
                              a general convention for the use of words:
                              word[0] stores the length of the word */

   matrix_TYP *erg,
             **ele,
              *tmp;

   erg = init_mat(G->gen_no * G->dim,1,"");
   for (i=0;i<G->gen_no;i++)
      Check_mat(map[i]);

   /* there are to cases: the function is called the first time
      for this generating set, then we have to calculate orbits and so on,
      otherwise we already got the desired words in words */
   if (word_flag == TRUE){
      found = (int *) calloc(G->gen_no+1 , sizeof(int));
      found++;
      ele = (matrix_TYP **) malloc(MIN_SPEICHER * sizeof(matrix_TYP *));
      ele_words = (int **) malloc(MIN_SPEICHER * sizeof(int*));

      for (i=0;i<G->gen_no;i++){
         ele[i] = map[i];
         found[red_pos(ele[i],G->gen,G->gen_no)] = TRUE;
         Check_mat(G->gen[i]);
         ele_words[i] = (int *) malloc(2 * sizeof(int));
         ele_words[i][0] = 1;
         ele_words[i][1] = i;
      }

      for (i=0;i<length && !and(found,G->gen_no);i++){

         for (j=0;j<number;j++){
            tmp = mat_mul(ele[i],map[j]);
            if (red_pos(tmp,ele,length) == (-1)){

               if (length >= speicher){
                  speicher = speicher + MIN_SPEICHER;
                  ele = (matrix_TYP **) realloc(ele,
                        speicher*sizeof(matrix_TYP *));
                  ele_words = (int **) realloc(ele_words,
                        speicher*sizeof(int *));
               }

               /* we found an new element of the orbit,
               set the matrix and the word */
               ele[length] = tmp;
               ele_words[length] = (int *) malloc((ele_words[i][0]+2)
                                               * sizeof(int*));
               memcpy(ele_words[length],ele_words[i],
                      (ele_words[i][0]+1) * sizeof(int));
               ele_words[length][0] = ele_words[i][0]+1;
               ele_words[length][ele_words[length][0]] = j;
               length++;

               /* it might be one of the generators we are looking for */
               found[red_pos(tmp,G->gen,G->gen_no)] = TRUE;
            }
            else{
               free_mat(tmp);
            }
         }
      }

      if (!and(found,G->gen_no)){
         fprintf(stderr,"reget_gen: shouldn't happen\n");
         exit(3);
      }

      /* now calculate the resulting cozycle alltogether */
      for (i=0;i<G->gen_no;i++){
         k = red_pos(G->gen[i],ele,length);
         tmp = ele[k];
         words[i] = ele_words[k];
         ele_words[k] = NULL;
         if (INFO_LEVEL & 16){
            printf("length %d \n",length);
            printf("red_pos %d \n",red_pos(G->gen[i],ele,length));
            Check_mat(tmp);
            put_mat(tmp,NULL,"tmp",2);
            put_mat(G->gen[i],NULL,"G->gen[i]",2);
            printf("word[%d]\n c ",i);
            for (j=1;j<=words[i][0];j++)
               printf("%d ",words[i][j]+1);
            printf("\n");
         }
         for (j=0;j<G->dim;j++){
            erg->array.SZ[i*G->dim +j][0] = tmp->array.SZ[j][G->dim];
         }
      }

      /* cleaning up memory */
      found--;
      free(found);
      for (i=G->gen_no;i<length;i++) free_mat(ele[i]);
      free(ele);
      for (i=0;i<length;i++) if (ele_words[i] != NULL) free(ele_words[i]);
      free(ele_words);
   }
   else{
      /* we are in the briliant position to know how the element looks
         like */
      for (i=0;i<G->gen_no;i++){
         tmp = copy_mat(map[words[i][1]]);
         for (j=2;j<=words[i][0];j++)
            mat_muleq(tmp,map[words[i][j]]);

         /* calculate the part of the cozycle belonging to the i-th
            generator */
         for (j=0;j<G->dim;j++)
            erg->array.SZ[i*G->dim +j][0] = tmp->array.SZ[j][G->dim];

         if (INFO_LEVEL & 16){
            printf("reget of %d-th generator\n",i+1);
            put_mat(tmp,NULL,NULL,2);
         }

         /* clean up tmp */
         free_mat(tmp);
      }
   }

   return erg;
} /* reget_gen(....) */

