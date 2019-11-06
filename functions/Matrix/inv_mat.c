#include "typedef.h"
#include "tools.h"
#include "matrix.h"
#include "longtools.h"


/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: inv_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
/*{{{}}}*/
/*{{{  rmat_inv       */
/*
*  result = rmat_inv( mat );
*
*  berechnet Inverse zu eine Nennermatrix. Das "r" steht fuer rational
*
*  matrix_TYP *result: Inverse zu mat, wird erzeugt
*  matrix_TYP *mat:    wird nicht veraendert
 */
static matrix_TYP *
rmat_inv (matrix_TYP *mat)                           
{  
rational **M, **II, *v, f ;
int n, i, j, k , l;
int *nuin_M, *nuin_I;
int h;
matrix_TYP *inv;

  f = Zero;
  n = mat->cols;
  nuin_M = (int *)calloc((n+1),sizeof(int));
  nuin_I = (int *)calloc((n+1),sizeof(int));
  
  /*
   * create rational version II of identity matrix        
   */
  /*{{{  */
  II = ( rational ** )malloc2dim( n, n, sizeof(rational) );
  memset2dim( (char ** )II, n, n, sizeof(rational), (char *)&Zero );
  for(i = 0; i < n; i++) {
    II[i][i] = One;
  }
  /*}}}  */
  /*
   * create rational version M of matrix mat
   */
  /*{{{  */
  M  = (rational **)malloc2dim( n, n, sizeof(rational) );
  for ( i  = 0; i  < n; i ++) {
    for ( j  = 0; j  < n; j ++) {
      M[i][j].z = mat->array.SZ[i][j];
      M[i][j].n = mat->array.N[i][j];
    }
  }
  /*}}}  */
  /*
   * Invert M by paralell transformation of M and II s.t M is unit 
   */                     

  /*
   *  first create triangular matrix
   */
  /*{{{  */
  for (i = 0; i < n; i++) {
    if ( M[i][i].z == 0 ) {
      /*
       * Swap a non-zero element into [i][i]-position        
       */
      /*{{{  */
      for (j = i+1; j < n && M[j][i].z == 0; j++ );
      if ( j == n ) {
        fprintf (stderr, "matrix is singular. Can't invert\n");
        exit (3);
      }
      v   =  M[i];  M[i] =  M[j];  M[j] = v;
      v   = II[i]; II[i] = II[j]; II[j] = v;
      /*}}}  */
    }
    /*
     * Find the non-zero entries in i-th row in M and II
     * 17.02.92, stored in var "f"
     */
    /*{{{  */
    nuin_M[0] =
    nuin_I[0] = 1;
    for(k = 0; k < n; k++) {
      if ( M[i][k].z != 0) nuin_M[nuin_M[0]++] = k;
      if (II[i][k].z != 0) nuin_I[nuin_I[0]++] = k;
    }
    f.z= M[i][i].z;
    f.n= M[i][i].n;
    /*}}}  */
    if ( (f.z != 1) || (f.n != 1) ) {
      /*
       * Normalize i-th row s.t. [i][i]-element is 1
       */
      /*{{{  */
      for (l = 1; l < nuin_M[0]; ++l) {
        k = nuin_M[l];
        M[i][k].z *= f.n;
        M[i][k].n *= f.z;
        Normal (&M[i][k]);
      }
      for (l = 1; l < nuin_I[0]; ++l) {
        k = nuin_I[l];
        II[i][k].z *= f.n;
        II[i][k].n *= f.z;
        Normal (&II[i][k]);
      }                              
      /*}}}  */
    }
    /*
     * Clear i-th column downwards
     */
    /*{{{  */
    for ( j = i+1; j < n; j++) {
      if ( M[j][i].z != 0 ) {
        f.z = M[j][i].z;
        f.n = M[j][i].n;
        for (l = 1; l < nuin_M[0]; ++l) {
          k = nuin_M[l];
          if(M[j][k].z != 0) {
            h= f.n * M[i][k].n;
            M[j][k].z *= h;
            M[j][k].z -= M[j][k].n * M[i][k].z * f.z;
            M[j][k].n *= h;
          } else {
            M[j][k].z= - M[i][k].z * f.z;
            M[j][k].n = M[i][k].n * f.n;
          }
          Normal (&M[j][k]);
        }
        for (l = 1; l < nuin_I[0]; ++l) {
          k = nuin_I[l];
          if(II[j][k].z != 0) {
            h= f.n * II[i][k].n;
            II[j][k].z *= h;
            II[j][k].z -= II[j][k].n * II[i][k].z * f.z;
            II[j][k].n *= h;
          } else {
            II[j][k].z = - II[i][k].z*f.z;
            II[j][k].n = II[i][k].n * f.n;
          }
          Normal (&II[j][k]);
        }
      }
    }
    /*}}}  */
  } /* M is now triangular */
  /*}}}  */
  /*
   * Clear columns upwards  
   */
  /*{{{  */
  for (i = n-1; i >= 0; i --) {
    nuin_I[0] = 1;
    for ( k = 0; k < n; k++) {
      if (II[i][k].z != 0) nuin_I[nuin_I[0]++] = k;
    }
    for ( j = i-1; j >= 0; j--) {
      if ( M[j][i].z != 0 ) {
        for (l = 1; l < nuin_I[0]; ++l) {
          k = nuin_I[l];
          h= M[j][i].n* II[i][k].n;
          II[j][k].z *= h;
          II[j][k].z -= M[j][i].z * II[j][k].n * II[i][k].z;
          II[j][k].n *= h;
          Normal (&II[j][k]);
        }                        
      }
    }
  }
  /*}}}  */
  /*
   * copy M to result matrix "inv"
   *
   */
  inv = init_mat(n,n,"r");
  for ( i  = 0; i < n  ; i++) {
    for ( j = 0; j < n; j++) {
      inv->array.SZ[i][j] = II[i][j].z;
      inv->array.N [i][j] = II[i][j].n;
    }
  }
  free2dim( (char **)M, n );
  free2dim( (char **)II, n );
  free(nuin_M);
  free(nuin_I);
  
  Check_mat( inv );

  return ( inv );
}

/*}}}  */
/*{{{  imat_inv*/
static matrix_TYP *
imat_inv (matrix_TYP *mat)
{  
matrix_TYP *inv, *help;

  /* don't relie on kgv2rat() and rat2kgv() being the inverse of each other
   * there might be integer overflows and we don't want to demolish the 
   * original matrix
   */
  help = copy_mat( mat );
  /* 
   * Allocate a real rational denominator matrix
   */
  kgv2rat( help );
  /*
   * compute the invers of it
   */
  inv = rmat_inv( help );    
  /*
   *  free workspace
   */
  free_mat( help );
  /*
   *  transform the result to lcm-representation. This might be stupid. 
   *  There could be integer overflow and in that case we'd better stay to the
   *  denominator matrix
   */
  rat2kgv( inv );       
  /*
   *  eigentlich sollte ein Check_mat() noetig sein. rat2kgv() berechnet 
   *  schon eine maximal gekuerzte Darstellung
   */                                         
#if 0
  Check_mat( inv );
#endif
  return ( inv );
}


/*{{{  pmat_inv, unchanged*/


/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *pmat_inv (mat)
@ matrix_TYP *mat;
@
@ calculates the inverse of mat modulo mat->prime
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
pmat_inv (matrix_TYP *mat)
{
int n, f, *v;
int i, j, k, **M, **II;
matrix_TYP *inv;

  n = mat->cols;
  init_prime(mat->prime);

  M = (int **)malloc2dim( n, n, sizeof(int) );
  memcpy2dim( (char **)M, (const char **)mat->array.SZ, n, n, sizeof(int) );
  II = (int **)calloc2dim( n, n, sizeof(int) );
  for ( i= 0; i< n; i ++) {
    II[i][i] = 1;
  }
  /*
   * Invert M by paralell transformation of M and II s.t M is unit 
   */
  for ( i= 0; i< n; i ++ ) {
    if ( M[i][i] == 0 ) {
   /*
    * Swap a non-zero element into [i][i]-position 
    */
      for ( j= i+1; j< n && M[j][i] == 0; j++ );
      if ( j == n ) {
        fprintf (stderr, "matrix is singular. Can't invert\n");
        exit (3);
      }
      v= M[i];   M[i] =  M[j];  M[j] = v;
      v= II[i]; II[i] = II[j]; II[j] = v;
    }
    f = M[i][i];
    if (f != 1) {
    /*
     * Normalize i-th row s.t. [i][i]-element is 1
     */
      for ( k= 0; k< n; k ++ ) {
        if ( k >= i ) {
          M[i][k] = P(M[i][k], - f);
        }
        II[i][k] = P(II[i][k], - f);
      }          
    }
    /*
     * Clear i-th column                     
     */
    for ( j= 0; j< n; j ++ ) {
      if ((j != i) && (f = M[j][i]) ) {
        for ( k= 0; k< n; k ++ ) {
          if (k > i) {
            M[j][k] = S(M[j][k], - P(f,M[i][k]));
          }
          II[j][k] = S(II[j][k], - P(f,II[i][k]));
        }
      }
    }
  }
  inv = init_mat(n,n,"p");
  free2dim( (char **)M, n );
  free2dim( (char **)inv->array.SZ, n  );
  inv->array.SZ = II;
  
  inv->prime = mat->prime;

  inv->flags.Symmetric = mat->flags.Symmetric;
  inv->flags.Diagonal  = mat->flags.Diagonal ;
  inv->flags.Scalar    = mat->flags.Scalar ;
  
  return (inv);
}

/*}}}  */
/*{{{  mat_inv*/
/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *mat_inv (mat)
@ matrix_TYP *mat;
@
@ calculates the inverse of mat
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
mat_inv (matrix_TYP *mat)
{
matrix_TYP *inv;

  if ( mat->rows != mat->cols ) {
    fprintf (stderr, "Can't invert non-square matrix\n");
    exit (3);
  }
  if ( mat->prime > 0 ) {
    inv = pmat_inv( mat );
  } else if ( mat->array.N != NULL ) {
    inv = rmat_inv( mat );
  } else {
    inv = imat_inv( mat );
  }
  return(inv);
}
/*}}}  */
