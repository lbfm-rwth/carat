#include"typedef.h"
#include"tools.h"
#include"matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: add_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *imat_add (L_mat, R_mat, Lc, Rc)
@ matrix_TYP *L_mat, *R_mat ;
@ int Lc, Rc;
@ 
@  calculates Lc * L_mat + Rc * R_mat
@---------------------------------------------------------------------------
@
\**************************************************************************/
/*{{{}}}*/
/*{{{  imat_add, exported*/
matrix_TYP *
imat_add (matrix_TYP *L_mat, matrix_TYP *R_mat, int Lc, int Rc)
{
matrix_TYP *S_mat;
int **L, **R, **S, rS, cS;
int i, j;

  rS = L_mat->rows;
  cS = L_mat->cols;
  S_mat = init_mat(rS,cS,"");
  
  S_mat->flags.Symmetric=(R_mat->flags.Symmetric&&L_mat->flags.Symmetric);
  S_mat->flags.Diagonal =(R_mat->flags.Diagonal && L_mat->flags.Diagonal);
  S_mat->flags.Scalar   =(R_mat->flags.Scalar   && L_mat->flags.Scalar);
  S_mat->flags.Integral =(R_mat->flags.Integral && L_mat->flags.Integral);

  /*
   *  ovfl_mul() is an expensive routine defined in ../tools/ovfl_mul.c
   *  It checks an integer-overflow using a division :(.
   *  But as the kgv should tend to grow fast, this seems to be saver
   */
  if(L_mat->kgv != 1) {
    Rc = ovfl_mul( Rc, L_mat->kgv );
  }
  if(R_mat->kgv != 1) {
    Lc = ovfl_mul( Lc, R_mat->kgv );
  }
  S_mat->kgv = ovfl_mul( L_mat->kgv , R_mat->kgv );
  
  S = S_mat->array.SZ;
  R = R_mat->array.SZ;
  L = L_mat->array.SZ;
#ifdef TEST_OVFL
  for ( i  = 0 ; i < rS; i++) {
    if(S_mat->flags.Diagonal) {
      S[i][i] = ovfl_mul( Lc , L[i][i] ) + ovfl_mul( Rc , R[i][i] );
    } else {
      for ( j  = 0 ; j  < cS  ; j ++  ) {
        S[i][j] = ovfl_mul( Lc , L[i][j] ) + ovfl_mul( Rc , R[i][j] );
      }
    }
  } 
#else
  for ( i  = 0 ; i < rS; i++) {
    if(S_mat->flags.Diagonal) {
      S[i][i] = Lc * L[i][i] + Rc * R[i][i];
    } else {
      for ( j  = 0 ; j  < cS  ; j ++  ) {
        S[i][j] = Lc * L[i][j] + Rc * R[i][j];
      }
    }
  }  
#endif
  Check_mat(S_mat);
  return(S_mat);
}

/*}}}  */
/*{{{  imat_addeq, exported*/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *imat_addeq (L_mat, R_mat, Lc, Rc)
@ matrix_TYP *L_mat, *R_mat ;
@ int Lc, Rc;
@ 
@  calculates L_mat = Lc*L_mat + Rc*R_mat 
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
imat_addeq (matrix_TYP *L_mat, matrix_TYP *R_mat, int Lc, int Rc)
{

int **L, **R, rS, cS;
int i, j;

  rS = L_mat->rows;
  cS = L_mat->cols;

  L_mat->flags.Symmetric=(R_mat->flags.Symmetric&&L_mat->flags.Symmetric);
  L_mat->flags.Diagonal =(R_mat->flags.Diagonal && L_mat->flags.Diagonal);
  L_mat->flags.Scalar   =(R_mat->flags.Scalar   && L_mat->flags.Scalar);
  L_mat->flags.Integral =(R_mat->flags.Integral && L_mat->flags.Integral);

  /*
   *  ovfl_mul() is an expensive routine defined in ../tools/ovfl_mul.c
   *  It checks an integer-overflow using a division :(.
   *  But as the kgv should tend to grow fast, this seems to be saver
   */
  if(L_mat->kgv != 1) {
    Rc = ovfl_mul( Rc, L_mat->kgv );
  }
  if(R_mat->kgv != 1) {
    Lc = ovfl_mul( Lc, R_mat->kgv );
  }
  L_mat->kgv = ovfl_mul( L_mat->kgv , R_mat->kgv );
  
  R = R_mat->array.SZ;
  L = L_mat->array.SZ;
#ifdef TEST_OVFL
  for ( i  = 0 ; i < rS; i++) {
    if(L_mat->flags.Diagonal) {
      L[i][i] = ovfl_mul( Lc , L[i][i] ) + ovfl_mul( Rc , R[i][i] );
    } else {
      for ( j  = 0 ; j  < cS  ; j ++  ) {
        L[i][j] = ovfl_mul( Lc , L[i][j] ) + ovfl_mul( Rc , R[i][j] );
      }
    }
  } 
#else
  for ( i  = 0 ; i < rS; i++) {
    if(L_mat->flags.Diagonal) {
      L[i][i] = Lc * L[i][i] + Rc * R[i][i];
    } else {
      for ( j  = 0 ; j  < cS  ; j ++  ) {
        L[i][j] = Lc * L[i][j] + Rc * R[i][j];
      }
    }
  }  
#endif
  Check_mat(L_mat);
  return(L_mat);
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *rmat_add (L_mat, R_mat, L_coeff, R_coeff)
@ matrix_TYP *L_mat, *R_mat ;
@ rational L_coeff, R_coeff;
@
@  calculates L_coeff * L_mat + R_coeff * R_mat
@---------------------------------------------------------------------------
@
\**************************************************************************/
/*}}}  */
/*{{{  rmat_add, exported*/
matrix_TYP *
rmat_add (matrix_TYP *L_mat, matrix_TYP *R_mat, rational L_coeff, rational R_coeff)
{
matrix_TYP *S_mat;
int **S, rS, cS;
int i, j;
rational Lc, Rc;

  Lc.z = L_coeff.z;
  Lc.n = L_coeff.n;
  Rc.z = R_coeff.z;
  Rc.n = R_coeff.n;
  Normal( &Lc );
  Normal( &Rc );

  rS = L_mat->rows;
  cS = L_mat->cols;
  S_mat = init_mat(rS,cS,"il");
  
  S_mat->flags.Symmetric=(R_mat->flags.Symmetric && L_mat->flags.Symmetric);
  S_mat->flags.Diagonal =(R_mat->flags.Diagonal && L_mat->flags.Diagonal);
  S_mat->flags.Scalar   =(R_mat->flags.Scalar   && L_mat->flags.Scalar);
  S_mat->flags.Integral =(R_mat->flags.Integral && L_mat->flags.Integral);
  if(!(S_mat->flags.Integral)) {
    int temp1= GGT(Lc.z, L_mat->kgv);
    int temp2= GGT(Rc.z, R_mat->kgv);
    Lc.z = Lc.z / temp1 * R_mat->kgv / temp2 * Rc.n;
    Normal(&Lc);
    Rc.z = Rc.z / temp2 * L_mat->kgv / temp1 * Lc.n;
    Normal(&Rc);
  
    S_mat->kgv = L_mat->kgv * R_mat->kgv * Lc.n * Rc.n;
    S_mat->kgv /= temp1;
    S_mat->kgv /= temp2;
    if (S_mat->kgv == 1) {
      S_mat->flags.Integral = TRUE;
    }
  }
  else {
    if((Lc.n == 1) && (Rc.n == 1)) {
      S_mat->kgv = 1;
    } else {
      S_mat->kgv = Lc.n * Rc.n ;
      S_mat->flags.Integral = FALSE;
      Lc.z *= Rc.n;
      Rc.z *= Lc.n;
    }
  }
  S= S_mat->array.SZ;
  for ( i  = 0 ; i < rS; i++) {
    if(S_mat->flags.Diagonal) {
      S[i][i] = Lc.z * L_mat->array.SZ[i][i] + Rc.z * R_mat->array.SZ[i][i];
    } else {
      for ( j  = 0 ; j  < cS  ; j ++  ) {
        S[i][j] = Lc.z * L_mat->array.SZ[i][j] + Rc.z * R_mat->array.SZ[i][j];
      }
    }
  }

  Check_mat(S_mat);
  return(S_mat);
}

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *rmat_addeq(L_mat, R_mat, L_coeff, R_coeff)
@ matrix_TYP *L_mat, *R_mat ;
@ rational L_coeff, R_coeff;
@
@  calculates L_coeff = L_coeff * L_mat + R_coeff * R_mat
@---------------------------------------------------------------------------
@
\**************************************************************************/
/*}}}  */
/*{{{  rmat_addeq, exported*/
matrix_TYP *
rmat_addeq (matrix_TYP *L_mat, matrix_TYP *R_mat, rational L_coeff, rational R_coeff)
{
int **L, **R, rS, cS;
int i, j;
rational Rc, Lc;

  rS = L_mat->rows;
  cS = L_mat->cols;
  Lc.z= L_coeff.z;
  Lc.n= L_coeff.n;
  Rc.z= R_coeff.z;
  Rc.n= R_coeff.n;
  Normal( &Lc );
  Normal( &Rc );

  L_mat->flags.Symmetric=(R_mat->flags.Symmetric && L_mat->flags.Symmetric);
  L_mat->flags.Diagonal =(R_mat->flags.Diagonal && L_mat->flags.Diagonal);
  L_mat->flags.Scalar   =(R_mat->flags.Scalar   && L_mat->flags.Scalar);
  L_mat->flags.Integral =(R_mat->flags.Integral && L_mat->flags.Integral);
  if(!(L_mat->flags.Integral)) {
    Lc.n= L_coeff.n * L_mat->kgv; 
    Normal(&Lc);
    Rc.n= R_coeff.n * R_mat->kgv; 
    Normal(&Rc);
    int temp1= GGT( Lc.n, Rc.n);
    Lc.z= Lc.z * Rc.n / temp1;
    Rc.z= Rc.z * Lc.n / temp1;
    L_mat->kgv= Lc.n *  Rc.n / temp1;
    if (L_mat->kgv == 1) {
      L_mat->flags.Integral = TRUE;
    }
  } else {
    if((L_coeff.n == 1) && (R_coeff.n == 1)) {
      L_mat->kgv = 1;
    } else {
      int temp1= GGT(Lc.n,Rc.n);
      L_mat->kgv = Lc.n * Rc.n / temp1;
      L_mat->flags.Integral = FALSE;
      Lc.z = Lc.z * Rc.n / temp1;
      Rc.z = Rc.z * Lc.n /temp1;
    }
  }
  L = L_mat->array.SZ;
  R = R_mat->array.SZ;
  for ( i  = 0 ; i < rS; i++) {
    if(L_mat->flags.Diagonal) {
      L[i][i] = L[i][i] * Lc.z + Rc.z * R[i][i];
    } else {
      for ( j  = 0 ; j  < cS  ; j ++  ) {
        L[i][j] = L[i][j] * Lc.z + Rc.z*R[i][j];
      }
    }
  }
  Check_mat(L_mat);
  return(L_mat);
}

/*}}}  */
/*{{{  pmat_add, exported*/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *pmat_add (L_mat, R_mat, L_coeff, R_coeff)
@ matrix_TYP *L_mat, *R_mat ;
@ int L_coeff, R_coeff;
@
@ calculates L_coeff * L_mat + R_coeff * R_mat modulo L_mat->prime
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
pmat_add (matrix_TYP *L_mat, matrix_TYP *R_mat, int L_coeff, int R_coeff)
{
matrix_TYP *S_mat;
int **L, **R, **A, rS, cS ;
int i, j;

  rS = L_mat->rows;
  cS = L_mat->cols;
  S_mat = init_mat(rS,cS,"p");
  S_mat->prime = L_mat->prime;
  S_mat->flags.Symmetric = (R_mat->flags.Symmetric&&L_mat->flags.Symmetric);
  S_mat->flags.Diagonal  = (R_mat->flags.Diagonal && L_mat->flags.Diagonal);
  S_mat->flags.Scalar    = (R_mat->flags.Scalar   && L_mat->flags.Scalar);
  A= S_mat->array.SZ;
  L= L_mat->array.SZ;
  R= R_mat->array.SZ;
  for ( i  = 0 ; i < rS; i++) {
    if(S_mat->flags.Diagonal) {
      A[i][i] = S(P(L_coeff,L[i][i]), P(R_coeff,R[i][i]));
    } else {
      for ( j  = 0 ; j  < cS  ; j ++  ) {
        A[i][j] = S(P(L_coeff,L[i][j]), P(R_coeff,R[i][j]));
      }
    }
  }
  Check_mat (S_mat);
  return(S_mat);
}

/*}}}  */
/*{{{  pmat_addeq, exported*/
/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *pmat_addeq (L_mat, R_mat, L_coeff, R_coeff)
@ matrix_TYP *L_mat, *R_mat ;
@ int L_coeff, R_coeff;
@
@ calculates L_mat = L_coeff * L_mat + R_coeff * R_mat modulo L_mat->prime
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
pmat_addeq (matrix_TYP *L_mat, matrix_TYP *R_mat, int L_coeff, int R_coeff)
{
int **L, **R, rS, cS ;
int i, j;

  rS = L_mat->rows;
  cS = L_mat->cols;
  L_mat->flags.Symmetric=(R_mat->flags.Symmetric&&L_mat->flags.Symmetric);
  L_mat->flags.Diagonal =(R_mat->flags.Diagonal && L_mat->flags.Diagonal);
  L_mat->flags.Scalar   =(R_mat->flags.Scalar   && L_mat->flags.Scalar);
  L= L_mat->array.SZ;
  R= R_mat->array.SZ;
  for ( i  = 0 ; i < rS; i++) {
    if(L_mat->flags.Diagonal) {
      L[i][i] = S(P(L_coeff,L[i][i]), P(R_coeff,R[i][i]));
    } else {
      for ( j  = 0 ; j  < cS  ; j ++  ) {
        L[i][j] = S(P(L_coeff,L[i][j]), P(R_coeff,R[i][j]));
      }
    }
  }
  Check_mat (L_mat);
  return(L_mat);
}

/*}}}  */
/*{{{  mat_add, exported*/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *mat_add (L_mat, R_mat, L_coeff, R_coeff)
@ matrix_TYP *L_mat, *R_mat ;
@ rational L_coeff, R_coeff;
@
@ calculates L_coeff * L_mat + R_coeff * R_mat
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
mat_add (matrix_TYP *L_mat, matrix_TYP *R_mat, rational L_coeff, rational R_coeff)
{
int Lc, Rc;
int lp;

  Lc = Rc = 1;
  if (( L_mat->cols != R_mat->cols) || (L_mat->rows != R_mat->rows )) {
    fprintf (stderr, "Can't add %dx%d with %dx%d\n",
                      L_mat->rows, L_mat->cols, R_mat->rows, R_mat->cols );
    exit (3);
  }
  if( L_mat->prime != 0 ) {
    if (L_mat->prime != R_mat->prime) {
      fprintf(stderr, "Error: Different primes in mat_add.\n");
      exit (3);
    }
    /*============================================================*\
    ||                                                            ||
    || The prime-field case:                                      ||
    || Initialize actuell prime and make sure that the            ||
    || Coeffizients are element of Fp.                            ||
    ||                                                            ||
    \*============================================================*/
    init_prime(L_mat->prime);
    lp = act_prime;
    if(L_coeff.n != 1) Lc = P(1,-L_coeff.n%lp);
    if(L_coeff.z != 1) Lc = P(Lc,L_coeff.z%lp);
    if(R_coeff.n != 1) Rc = P(1,-R_coeff.n%lp);
    if(R_coeff.z != 1) Rc = P(Rc,R_coeff.z%lp);
    return  pmat_add(L_mat,R_mat,Lc,Rc);
  } else {
    if ( L_mat->array.N != NULL ) {
      rat2kgv( L_mat );
    }    
    if ( R_mat->array.N != NULL ) {
      rat2kgv( R_mat );
    }
    return rmat_add(L_mat,R_mat,L_coeff, R_coeff);
  }
}

/*}}}  */
/*{{{  mat_addeq, exported*/
/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *mat_addeq (L_mat, R_mat, L_coeff, R_coeff)
@ matrix_TYP *L_mat, *R_mat ;
@ rational L_coeff, R_coeff;
@
@ calculates L_mat = L_coeff * L_mat + R_coeff * R_mat
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
mat_addeq (matrix_TYP *L_mat, matrix_TYP *R_mat, rational L_coeff, rational R_coeff)
{
int Lc, Rc;
int lp;

  Lc = Rc = 1;
  if (( L_mat->cols != R_mat->cols) || (L_mat->rows != R_mat->rows )) {
    fprintf (stderr, "Can't add %dx%d with %dx%d\n",
                      L_mat->rows, L_mat->cols, R_mat->rows, R_mat->cols );
    exit (3);
  }
  if( L_mat->prime != 0 ) {
    if (L_mat->prime != R_mat->prime) {
      fprintf(stderr, "Error: Different primes in mat_add.\n");
      exit (3);
    }
    /*============================================================*\
    ||                                                            ||
    || The prime-field case:                                      ||
    || Initialize actuell prime and make sure that the            ||
    || Coeffizients are element of Fp.                            ||
    ||                                                            ||
    \*============================================================*/
    init_prime(L_mat->prime);
    lp = act_prime;
    if(L_coeff.n != 1) Lc = P(1,-L_coeff.n%lp);
    if(L_coeff.z != 1) Lc = P(Lc,L_coeff.z%lp);
    if(R_coeff.n != 1) Rc = P(1,-R_coeff.n%lp);
    if(R_coeff.z != 1) Rc = P(Rc,R_coeff.z%lp);
    return pmat_addeq(L_mat,R_mat,Lc,Rc);
  } else {
    if ( L_mat->array.N != NULL ) {
      rat2kgv( L_mat );
    }    
    if ( R_mat->array.N != NULL ) {
      rat2kgv( R_mat );
    }
    return rmat_addeq(L_mat,R_mat,L_coeff, R_coeff);
  }
}

/*}}}  */
