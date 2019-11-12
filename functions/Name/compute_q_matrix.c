#include "typedef.h"
#include "utils.h"
#include "matrix.h"
#include "idem.h"
#include "bravais.h"
#include "orbit.h"
#include "sort.h"

/*Wahrscheinlich verursacht das Rechnen mit Matrizen, die nicht Ganzzahlig sind
  Probleme. */

/* This function compares two matrices in lexicographic order.
   It is supposed that they behave well (cols etc equal). */
static int cmp_mat_lex (matrix_TYP *A, matrix_TYP *B)
{
  int i, j;

  
      for (i=0; i<A->rows; i++) {
	for (j=0; j<A->cols; j++) {
	  if (A->array.SZ[i][j] != B->array.SZ[i][j]) {
	    if (A->array.SZ[i][j] < B->array.SZ[i][j])
	      return -1;
	    else
	      return  1;
	  }
	}
      }
      
      return 0;
}

static void sort_cols (matrix_TYP *N)
{
  
  int i, j, is_different, swap_cols;
  
  for (i=1; i<N->cols; i++)
    {
      is_different = FALSE, swap_cols = FALSE;
      
      for (j=0; j<N->rows && is_different == FALSE; j++)
	if (N->array.SZ[j][i] > N->array.SZ[j][i-1])
	  is_different = TRUE;
	else if (N->array.SZ[j][i] < N->array.SZ[j][i-1])
	  is_different = TRUE,
	    swap_cols = TRUE;
      
      if (swap_cols)
	{
	  col_per (N, i-1, i);
	  i = 0;
	}
      
    }
  
}

static void permute_and_sort (matrix_TYP *N, int number)
{
   int i, j, k = 1;
   matrix_TYP **set;
   int used = -1, unused = 0;
   
if (number != 1)
{
   for (i=2; i<=number; i++)
     k *= i;

   set = (matrix_TYP **)xmalloc((k+1) * sizeof(matrix_TYP *));

   set [0] = copy_mat (N);

   while (used < unused)
     {

	used ++;

       set [unused + 1] = copy_mat (set [used]);
       
       row_per (set [unused + 1], N->rows - number, N->rows - number + 1);
       
       for (i=0; i<=unused && cmp_mat_lex (set [unused + 1], set [i]) != 0; i++)
	 ;
       
       
       if (i == unused + 1)
	 unused ++;
       else
	 free_mat (set [unused + 1]);
       
       
       
       set [unused + 1] = copy_mat (set [used]);
       
       row_per (set [unused + 1], N->rows - number, N->rows - 1);
       
       for (i=1; i <= number - 1; i++)
	 row_per (set [unused + 1], N->rows - i - 1, N->rows - i);
       
       for (i=0; i<=unused && cmp_mat_lex (set [unused + 1], set [i]) != 0; i++)
	 ;
       
       if (i == unused + 1)
	 unused ++;
       else
	 free_mat (set [unused + 1]);
       
     }
   
   
   for (i=0; i<=used; i++)
     sort_cols (set[i]);
   
   
   k = 0;
   for (i=1; i<=used; i++)
     if (cmp_mat_lex (set [k], set [i]) > 0)
       k = i;
   
   
   for (i=0; i<N->rows; i++)
     for (j=0; j<N->cols; j++)
       N->array.SZ [i][j] = set [k]->array.SZ [i][j];
   
   for (i=0; i<=unused; i++)
     free_mat (set [i]);
   
   free (set);
}
else 
{
sort_cols (N);
}
}  

static int order(matrix_TYP *A)
{
  
  int i = 1;
  
  matrix_TYP *B;
  
  B = copy_mat(A);
  Check_mat(B);
  
  while (! ( B->flags.Scalar && B->array.SZ[0][0] == 1)){
    mat_muleq(B,A);
    Check_mat(B);
    i++;
  }
  
  free_mat(B);
  return i;
}

matrix_TYP *compute_q_matrix (bravais_TYP *G)
{
  matrix_TYP **Group, **Idem, **representant, *REP;
  matrix_TYP *I, *F, *output_mat;
  matrix_TYP *waste, *waste2;
  int orb_alg_options[6];
  int group_order, *length_conj_class;
  /* int reserved_conj_class = 0; */
  int class_counter = 0;
  int idem_no;
  int i, j, h;


  I = einheitsmatrix( G->dim );
  F = rform(G->gen, G->gen_no, I, 101);
  Idem = idempotente(G->gen, G->gen_no, F, &idem_no, &i, &j);
  free_mat(F);
  for (h=0; h<i+j; h++)
    free_mat(Idem[idem_no+h]);
  /* Here we free some Idem-elements that are not needed for our aims. */

  /* As a first step compute the whole group out of the
     Generators given by the "bravais_TYP *G".
     This is done by getting the orbit of the unit element
     under operation of the group.*/

    orb_alg_options[0] = 0;
    orb_alg_options[1] = 0;
    orb_alg_options[2] = 0,
    orb_alg_options[3] = 0,
    orb_alg_options[4] = 0,
    orb_alg_options[5] = 0;
  
  Group = orbit_alg( I, G, NULL, orb_alg_options, &group_order);
  mat_quicksort(Group, 0 , group_order - 1,mat_comp);

  /* The next step simply consists in dividing the whole group into its
     conjugacy classes. This is done by getting the orbit under conjugation
     for every group element that has not yet been found to lie in another
     conjugacy class.
     Technically this is done by setting up a list "not_yet_found" that
     contains a flag for every group element saying if it was already found.
     Das Programm startet am Anfang der Liste der Gruppenelemente und berechnet
     nach einander von ...*/
  
  
    orb_alg_options[0] = 4;
    orb_alg_options[1] = 0;
    orb_alg_options[2] = 0;
    orb_alg_options[3] = FALSE;
    orb_alg_options[4] = 0;
    orb_alg_options[5] = 0;
  

  REP = orbit_representatives(Group,
                              group_order,
                              G,
                              orb_alg_options,
                              &class_counter,
                              1);

  representant = (matrix_TYP **)xmalloc(class_counter * sizeof(matrix_TYP *));

  for (i=0;i<class_counter;i++){
     representant[i] = Group[REP->array.SZ[0][i]];
     Group[REP->array.SZ[0][i]] = NULL;
  }
  length_conj_class = REP->array.SZ[1];
  REP->array.SZ[1] = (int *) malloc(1*sizeof(int));

  for(i=0; i<group_order; i++){
    if (Group[i] != NULL) free_mat(Group[i]);
  }
  free(Group);
  free_mat(REP);
 
  output_mat = init_mat( 9 + idem_no, class_counter, "");
  
  /* Achtung: Probleme bei nicht ganzzahligen Matrizen !!! */

  /* Put the length of the conjugation classes into the 0-th row of the
     output matrix. */
  for( i=0; i<class_counter; i++)
    output_mat->array.SZ[0][i] = length_conj_class[i];
  free(length_conj_class);
  
  /* Put the order of the representant into the 1-st row of the output matrix*/
  for( i=0; i<class_counter; i++)
    output_mat->array.SZ[1][i] = order( representant[i] );
  
  /* Put the trace of the conjugation classes into the 2-nd row of the
     output matrix. */
  for( i=0; i<class_counter; i++)
    output_mat->array.SZ[2][i] = trace( representant[i] );
  
  
  /* Put the trace of the square of elements of every conjugation classes
     into the 3-nd row of the output matrix. */
  for( i=0; i<class_counter; i++)
    {
      waste = mat_mul( representant[i], representant[i]);
      
      output_mat->array.SZ[3][i] = trace(waste);
      
      free_mat(waste);
    }
  
  /* Put the trace of the cube of elements of every conjugation classes
     into the 4-th row of the output matrix. */
  for( i=0; i<class_counter; i++)
    {
      waste2 = mat_mul( representant[i], representant[i]);
      waste = mat_mul( waste2, representant[i] );
      
      output_mat->array.SZ[4][i] = trace( waste );
      
      free_mat(waste);
      free_mat(waste2);
    }
  
  /* Get the rank of the fixspace of every element */
  for (i=0; i<class_counter; i++)
    {
      waste = copy_mat( representant[i] );
      
      for (j=0; j<I->rows; j++)
	waste->array.SZ[j][j]--; /* This means: matrix minus identity */ 
      
      Check_mat( waste );
      output_mat->array.SZ[5][i] = tgauss( waste );
      free_mat( waste );
    }
  
  /* Compute the rank of the antifixspace of every element. */
  for (i=0; i<class_counter; i++)
    {
      waste = copy_mat( representant[i] );
      
      for (j=0; j<I->rows; j++)
	waste->array.SZ[j][j]++;/* This means: matrix plus identity */
      
      Check_mat( waste );
      output_mat->array.SZ[6][i] = tgauss( waste );
      free_mat( waste );
    }

  /* get the rank of (x-1)^2 */
  for (i=0; i<class_counter; i++)
    {
      waste = copy_mat( representant[i] );
      for (j=0; j<I->rows; j++)
	waste->array.SZ[j][j]--;
      
      Check_mat( waste );
      mat_muleq(waste,waste);
      Check_mat( waste );
      output_mat->array.SZ[7][i] = tgauss( waste );
      free_mat( waste );
    }

  /* put the rank of (x+1)^2 into the 8-th row of the output_mat*/
   for (i=0; i<class_counter; i++)
     {
       waste = copy_mat( representant[i] );

       for (j=0; j<I->rows; j++)
	 waste->array.SZ[j][j]++;

       Check_mat( waste );
       mat_muleq( waste, waste );
       Check_mat( waste );
       output_mat->array.SZ[8][i] = tgauss( waste );
       free_mat( waste );
     }

   /* Write into the 9-th row of output_mat the trace of idem??? */
   for (i=0; i<class_counter; i++)
     for (j=0; j<idem_no; j++)
       {
	 waste = mat_mul(Idem[j],representant[i]);
	 output_mat->array.SZ[9+j][i] = trace(waste)/waste->kgv;
	 free_mat(waste);
       }
   for(i=0; i<idem_no; i++)
     free_mat(Idem[i]);
   free(Idem);
   
   /* The matrix has now to be ordered by shuffling columns and rows. */
   permute_and_sort (output_mat, idem_no);
   
   
   
  
   free_mat(I);
   
   for(i=0; i<class_counter; i++)
     free_mat(representant[i]);
   free(representant);   

  return output_mat;
}
