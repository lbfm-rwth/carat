#include "typedef.h"
#include "matrix.h"
#include "gmp.h"
#include "zass.h"
#include "getput.h"
#include "tools.h"
#include "contrib.h"

/**************************************************************************
@
@--------------------------------------------------------------------------
@
@ FILE: torsionfree.c
@
@--------------------------------------------------------------------------
@
***************************************************************************/

/* compare two matrix */
static int is_equal(matrix_TYP *x,
                    matrix_TYP *y)
{
  int i,
      j,
      cols,
      rows;

  cols = x->cols;
  rows = x->rows;

  for (i = 0 ; i < rows ; i++){
    for(j = 0 ; j < cols ; j++){
      if ((y->kgv * x->array.SZ[i][j]) != (x->kgv * y->array.SZ[i][j])){
	 return FALSE;
      }
    }
  }
       
  return TRUE;
}   


/* compare two matrix (linear part only) */
static int lin_is_equal(matrix_TYP *x,
                        matrix_TYP *y)
{
  int i,
      j,
      cols,
      rows;

  cols = x->cols - 1;
  rows = x->rows - 1;

  for (i = 0 ; i < rows ; i++)
    for(j = 0 ; j < cols ; j++)
       if ((x->array.SZ[i][j] * y->kgv) != (y->array.SZ[i][j] * x->kgv))
	 return FALSE;
  return TRUE;
}


/* calculate the linear part P of R */
matrix_TYP *lin(matrix_TYP *x)
{
  matrix_TYP *y;
  y = copy_mat (x);
  real_mat (y, y->rows-1, y->cols-1);
  Check_mat (y);
  return y;
}

/* calculate the cyclotomic polyn */
static matrix_TYP *pol_cycl(matrix_TYP *x,
                           int p)
{
  rational one;

  int i;

  matrix_TYP *y,
             *Id,
             *RES;

  one.z = 1;
  one.n = 1;

  y = copy_mat(x);
  Id = init_mat(x->rows, x->cols, "i1");
  RES = init_mat(x->rows, x->cols, "i1");

  for(i=1 ; i < p ; i++){
    RES = mat_muleq( RES , y );
    RES = mat_addeq( RES, Id , one , one);
  }
  free_mat (y);
  free_mat (Id);
  
  return RES;
}

/* calculate the n-power of a matrix */
static matrix_TYP *power(matrix_TYP *x,
                         int n)
{
  int i;

  matrix_TYP *y;

  y = init_mat(x->rows, x->cols, "i1");
  if (n == 0){
    return y;
  }
  else{
  for(i=0 ; i < n ; i++){
    y = mat_muleq( y , x );
  }
  
  return y;
  
  }
}


/* calculate the order of a matrix */
static int ORDER(matrix_TYP *A)
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


/**************************************************************************
@
@--------------------------------------------------------------------------
@
@ int *torsionfree(bravais_TYP *R,
@                  int *order_out,
@                  int *number_of_conjugacy_classes){
@
@  decides whether the space group in R is torsion free, and ONLY
@  IN THIS CASES decides whether the group has trivial center or
@  not.
@  The result is a pointer A, say, which points to two integers with
@  the following convention:
@  A[0] == TRUE iff the group R is torsion free.
@  A[1] == TRUE iff A[0] == TRUE and the group R has trivial center.
@
@  bravais_TYP *R: a bravais_TYP with generators which together
@                  with Z^n generate the space group in Question.
@  int *order    : The order of the point group is stored here after
@                  the calculation has been done
@  int *number_of_conjugacy_classes: The number of conjugacy classes of
@                                    the point group of R is stored here.
@--------------------------------------------------------------------------
@
***************************************************************************/
int *torsionfree(bravais_TYP *R,
                 int *order_out,
                 int *number_of_conjugacy_classes){

int i,
    j,
    k,
    l,
    new_index,
    act_number,
    dim,
   *mirror,
   *processed,
   *list_primes,
    num_gen,
    order,
   *RES;


matrix_TYP **ele,
           **gen,
           **repres,
           **gen_inv,
           **gen_ident,
            *conj,
            *new_elem,
            *M,
            *X,
            *linear,
            *w,
            *y,
            *bahn,
            *bahn1,
            *bahn2,
            *ed1,
            *ed2;

RES = (int *) calloc(2 , sizeof(int));

num_gen = R->gen_no;
gen = R->gen;
dim = R->dim;


ele = (matrix_TYP **) malloc(110000 * sizeof(matrix_TYP *));
gen_inv = (matrix_TYP **) malloc(num_gen * sizeof(matrix_TYP *));
repres = (matrix_TYP **) malloc(110000 * sizeof(matrix_TYP *));
gen_ident = (matrix_TYP **) malloc(num_gen * sizeof(matrix_TYP *));

/* We will generate a list ele[] of elements of R that are pre image P  */
/* "order" is the order of P */

ele[0] = init_mat(dim, dim, "i1");

order = 1;

 for(i=0; i < order ; i++){
   for(j=0 ; j < num_gen ; j++){
     new_elem = mat_mul(ele[i] , gen[j]);
     Check_mat(new_elem);
     for(k=0 ; k < order && lin_is_equal(new_elem, ele[k])==FALSE; k++);
     if (k == order){
       ele[order] = new_elem;
       order++;
     }
     else{
       free_mat(new_elem);
     }
   }
 }
 
/* If P is trivial, then it's torsion-free ( Z^n)  */
 
 if (order == 1){  
   RES[0] = TRUE;
   order_out[0] = 1;
   number_of_conjugacy_classes[0] = 1;
   free_mat (ele[0]);
   free (ele);
 }
 else{
   /* for (i=0 ; i< order ; i++)
     put_mat(ele[i] , NULL , "" , 2); */

     order_out[0] = order;

   /* Now we calculate the (pre-images of) conjugacy classes of elements of P */

   mirror = (int *) malloc(order * sizeof(int));
   processed = (int *) malloc(order * sizeof(int));


   act_number = 0 ;

   for( i=0 ; i < order ; i++){
     mirror[i] = -1;
     processed[i] = 0;
   }
 
   for (i=0; i < num_gen ; i++){
     gen_inv[i] = mat_inv(gen[i]);
   }
 
   for (i=0 ; i < order ; i++){
     if (mirror[i] == -1){   
       act_number++;
       mirror[i] = act_number;
     
       for ( j = i ; j < order ; j++){
	 new_index = j;
	 if (mirror[j] == act_number && processed[j] == 0){
	   processed[j] = 1;

	   for (k=0 ; k < num_gen ; k++){
	     conj = mat_kon( gen[k] , ele[j] , gen_inv[k]);
	     for(l=0 ; l < order && lin_is_equal(conj, ele[l])==FALSE; l++);
	     mirror[l] = act_number;
	     if (l < new_index && processed[l] == 0){
	       new_index = l - 1;
	     }
	     free_mat(conj);
	   }
	 }
	 j = new_index;
       }
     }
   }
   number_of_conjugacy_classes[0] = act_number;

   for (i=0; i<num_gen; i++)
     free_mat (gen_inv[i]);
   free (gen_inv);

   /* we will generate a list repres[] of repres of conjug classes */
   /* act_number is the number of classes. */

   if(act_number == order){
     for (i=0 ; i < order ; i++)
       repres[i] = ele[i];
   }
   else{
     for( i=0 ; i < act_number ; i++){
       for(j = 0 ; j < order && mirror[j] != i+1 ; j++);
       repres[i] = ele[j];
     }
   }
  
   free (mirror); free (processed);

   list_primes = factorize(order);

   for (i=2 ; i < 98 ; i++){
     if ( list_primes[i] != 0){
       for (j = 0 ; j < act_number ; j++){
	 linear = lin(repres[j]);
	 if (i == ORDER (linear)){
	   
	   bahn = pol_cycl(linear, i);
	   bahn1 = tr_pose (bahn);

	   free_mat (bahn);
           free_mat (linear);

	   w = power (repres[j],i);
	   y = init_mat(1,dim - 1,"");
	   
	   for (l=0; l < dim-1; l++){
	     y->array.SZ[0][l] = w->array.SZ[l][dim-1];
	   }
	  	   
	   bahn2 = copy_mat (bahn1);
	   real_mat (bahn2 , dim , dim-1); 
	   for (k=0; k < dim-1; k++){
	     bahn2->array.SZ[dim-1][k] = y->array.SZ[0][k];
	   }
	  
	   ed1 = elt_div(bahn1); 
	   ed2 = elt_div(bahn2);

	   if(is_equal(ed1 , ed2) == TRUE){
	     free (list_primes);
             free (gen_ident);
	     free_mat (bahn1); free_mat (ed1);
	     free_mat (bahn2); free_mat (ed2);
	     free_mat (y); free_mat (w);
	     for (i=0; i<order; i++)
	       free_mat (ele[i]);
	     free (ele);
	     free (repres);
	     return RES;
	   }
	   free_mat (bahn1); free_mat(ed1);
	   free_mat (bahn2); free_mat(ed2);
	   free_mat (y); free_mat (w);
	 }
         else{
    	   free_mat (linear);
         }
        }
     }
   }
   RES[0] = TRUE;



   for (i=0; i<order; i++)
     free_mat (ele[i]);
   free (ele);
   free (repres);
   free (list_primes);


 /**  CALCULATE THE CENTRE **/
 
  /* We look for the module L^P, the elements of the lattice
       L centralized  by P, the point-group   */

  /* We create a list of matrix, gen_ident, equal to linear part
              of gen minus the Identity  */ 

  for(i=0 ; i < num_gen ; i++){
    gen_ident[i] = init_mat(dim-1,dim-1,"");
    for (j=0 ; j < dim -1 ; j++){
      for (k = 0 ; k < dim - 1 ; k++){
	if (k == j ){
	  gen_ident[i]->array.SZ[j][k] = (gen[i]->array.SZ[j][k]/gen[i]->kgv) - 1;
	}
	else{
	  gen_ident[i]->array.SZ[j][k] = gen[i]->array.SZ[j][k]/gen[i]->kgv;
	}
      }
    }
  }
  

 /* construct a greater matrix M with the ones above  */
  
  
  M = init_mat (num_gen*(dim-1) , dim-1 , "");
  
  for (i = 0 ; i < num_gen ; i++){
    for (j = i*(dim-1) ; j < (i+1)*(dim-1) ; j++){
      for (k = 0 ; k < dim - 1 ; k++){
	M->array.SZ[j][k] = gen_ident[i]->array.SZ[j-i*(dim-1)][k];
      }
    }
  }
    
  for (i=0; i<num_gen; i++)
    free_mat (gen_ident[i]);
  free (gen_ident);
  
   X =  solve_mat(M);
    /* the rows of X are Z-basis for the solution of MX = 0  */
  free_mat (M);
  
  if (X->rows == 0){
     RES[1] = TRUE;
  }
  else if (X->rows == 1){
    for (i = 0 ; i < dim-1 && X->array.SZ[0][i] == 0; i++);
    if (i == dim-1){
      RES[1] = TRUE;
    }
  }
  free_mat (X);
  
 }

  return RES;
}
