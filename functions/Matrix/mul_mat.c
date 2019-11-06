#include "typedef.h"
#include "tools.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: mul_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
/*{{{}}}*/
/*{{{  rmat_mul, done*/
static matrix_TYP *
rmat_mul (matrix_TYP *L_mat, matrix_TYP *R_mat)

{
matrix_TYP *P_mat;
int **L, **LN, **R, **RN, **P, **PN ;
int  rL, cL, cR, rP, cP ;
int i, j, k ; 
int tmp_ggt, tmp_z, tmp_n, tmp_kgv;

  cL =        L_mat->cols;
  rL = rP =   L_mat->rows;
  cR = cP =   R_mat->cols;
  
  P_mat = init_mat(rP,cP,"r");
  
  L = L_mat->array.SZ;
  R = R_mat->array.SZ;
  P = P_mat->array.SZ;
  LN = L_mat->array.N;
  RN = R_mat->array.N;
  PN = P_mat->array.N;
  
  for (i = 0; i < rL; i++)
    for (j = 0; j < cL; j++ ) 
    {
      if ( L[i][j]) {
        if ( R_mat->flags.Diagonal ) {
          if (R[j][j]) { 
            P[i][j]  =  L[i][j] *  R[j][j];
            PN[i][j] = LN[i][j] * RN[j][j];
          }
        } else {
          for (k = 0; k < cR; k++) {
            if (R[j][k]) { 
              tmp_z =  L[i][j] *  R[j][k];
              tmp_n = LN[i][j] * RN[j][k];
              /*
               *   direkt kuerzen! Ueberlaufgefahr
               */
              Normal2( &tmp_z, &tmp_n );
              tmp_ggt = GGT( PN[i][k], tmp_n );
              if ( PN[i][k] > tmp_n ) {
                tmp_kgv = PN[i][k] / tmp_ggt * tmp_n;
              } else {
                tmp_kgv = tmp_n/ tmp_ggt * PN[i][k];
              }                                   
              P[i][k]  = P[i][k] * (tmp_kgv / PN[i][k]);
              P[i][k] += tmp_z * (tmp_kgv / tmp_n );
              PN[i][k] = tmp_kgv;
              Normal2( &P[i][k], &PN[i][k]);
            }
          }
        }           
      }
    }
  
  P_mat->kgv = 1;
  P_mat->flags.Integral  = FALSE;
  P_mat->flags.Scalar = L_mat->flags.Scalar && R_mat->flags.Scalar;
  P_mat->flags.Diagonal = L_mat->flags.Diagonal && R_mat->flags.Diagonal;
  P_mat->flags.Symmetric = P_mat->flags.Diagonal;
  Check_mat( P_mat);
  return ( P_mat );
}
/*}}}  */
/*{{{  intmat_mul, done*/
static matrix_TYP *
intmat_mul (matrix_TYP *L_mat, matrix_TYP *R_mat)

{
matrix_TYP *P_mat;
int **L, **R, **P ;
int *L_i,
    *R_j,
    *P_i;
int  rL, cL, cR, rP, cP ;
int i, j, k ;

  cL =        L_mat->cols;
  rL = rP =   L_mat->rows;
  cR = cP =   R_mat->cols;
  
  P_mat = init_mat(rP,cP,"");
  
  L = L_mat->array.SZ;
  R = R_mat->array.SZ;
  P = P_mat->array.SZ;
  
  for (i = 0; i < rL; i++){
    L_i = L[i];
    P_i = P[i];
    for (j = 0; j < cL; j++) 
    {
      R_j = R[j];
      if ( L_i[j]) {
        if ( R_mat->flags.Diagonal ) 
        {
          if (R_j[j]) P_i[j] = L_i[j] * R_j[j];
        }
        else
        {
          for (k = 0; k < cR; k++)
            if (R_j[k])
               P_i[k] += L_i[j] * R_j[k];
        }
      }
    }
  } 
  P_mat->kgv = L_mat->kgv * R_mat->kgv;
  P_mat->flags.Integral  = (P_mat->kgv == 1);
  P_mat->flags.Scalar = L_mat->flags.Scalar && R_mat->flags.Scalar;
  P_mat->flags.Diagonal = L_mat->flags.Diagonal && R_mat->flags.Diagonal;
  P_mat->flags.Symmetric = P_mat->flags.Diagonal;
  Check_mat( P_mat);
  return ( P_mat );
}
/*}}}  */
/*{{{  pmat_mul, done*/


/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *pmat_mul (L_mat, R_mat)
@ matrix_TYP *L_mat, *R_mat ;
@
@ calculates L_mat * R_mat modulo L_mat->prime
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
pmat_mul (matrix_TYP *L_mat, matrix_TYP *R_mat)
{
matrix_TYP *P_mat;

int **L, rL, cL, cR,  **R, **M ;

int i, j, k;
	
rL = L_mat->rows;
cL = L_mat->cols;
cR = R_mat->cols;
P_mat = init_mat(rL,cR,"p");
L = L_mat->array.SZ;
R = R_mat->array.SZ;
M = P_mat->array.SZ;
P_mat->prime = L_mat->prime;
init_prime(P_mat->prime);

for ( i = 0 ; i < rL; i++)
	for (j = 0; j < cL; j++ )
		if ( L[i][j])
			for (k = 0 ; k < cR; k++)
				if ( R[j][k] )
					M[i][k] = S(M[i][k],P(L[i][j],R[j][k]));
Check_mat(P_mat);
return ( P_mat );
}

/*}}}  */
/*{{{  mat_mul, done*/


/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *mat_mul (L_mat, R_mat)
@ matrix_TYP *L_mat, *R_mat ;
@
@ calculates L_mat * R_mat.
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
mat_mul (matrix_TYP *L_mat, matrix_TYP *R_mat)

{
matrix_TYP *P_mat, *tmp_mat;


  if ( L_mat->cols != R_mat->rows ) {
    fprintf (stderr, "Can't multiply %dx%d with %dx%d\n", L_mat->rows,
                                L_mat->cols, R_mat->rows, R_mat->cols );
    exit (3);
  }
  if ( L_mat->prime != R_mat->prime ) {
    fprintf(stderr,"Can`t multiply Fp-matrix and Z-matrix\n");
    exit(3);
  }
  if ( L_mat->prime ) {
    P_mat = pmat_mul(L_mat, R_mat);
  } else {
    if ( L_mat->array.N != NULL ) {
      if ( R_mat->array.N != NULL ) {
        P_mat = rmat_mul( L_mat, R_mat );
      } else {
        tmp_mat = copy_mat( R_mat );
        kgv2rat( tmp_mat );
        P_mat = rmat_mul( L_mat, tmp_mat );
        free_mat( tmp_mat );
      } 
    } else if ( R_mat->array.N != NULL ) {
      tmp_mat = copy_mat( L_mat );
      kgv2rat( tmp_mat );
      P_mat = rmat_mul( tmp_mat, R_mat );
      free_mat( tmp_mat );
    } else {
      P_mat = intmat_mul(L_mat, R_mat);
    }
  }
  return(P_mat);
}

/*}}}  */
/*{{{  mat_muleq, unchanged*/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *mat_muleq(L_mat,R_mat)
@ matrix_TYP *L_mat, *R_mat;
@
@ calculates L_mat = L_mat * R_mat.
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
mat_muleq (matrix_TYP *L_mat, matrix_TYP *R_mat)

{
matrix_TYP *H_mat;
matrix_TYP merk;

  H_mat = mat_mul(L_mat,R_mat);

  memcpy( &merk, L_mat, sizeof(matrix_TYP) );
  memcpy( L_mat, H_mat, sizeof(matrix_TYP) );

  /* inserted 21/1/97 tilman */
  memcpy( H_mat, &merk, sizeof(matrix_TYP) );


  free_mat(H_mat);

  return(L_mat);
}
/*}}}  */
/*{{{  mat_kon*/


/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *mat_kon(L_mat,M_mat,R_mat)
@ matrix_TYP *L_mat, *M_mat, *R_mat;
@
@ calculates L_mat * M_mat * R_mat
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *
mat_kon (matrix_TYP *L_mat, matrix_TYP *M_mat, matrix_TYP *R_mat)
{
matrix_TYP *P_mat, *T_mat;


  if ( L_mat->cols != M_mat->rows ) 
  {
    fprintf (stderr, "Can't multiply %dx%d with %dx%d\n", L_mat->rows,
                                L_mat->cols, M_mat->rows, M_mat->cols );
    exit (3);
  }
  if ( M_mat->cols != R_mat->rows ) 
  {
    fprintf (stderr, "Can't multiply %dx%d with %dx%d\n", L_mat->rows,
                                L_mat->cols, R_mat->rows, R_mat->cols );
    exit (3);
  }
  if ( L_mat->prime != R_mat->prime )
  {
    fprintf(stderr,"Can`t multiply Fp-matrix and Z-matrix\n");
    exit(3);
  }
  if ( L_mat->prime != M_mat->prime )
  {
    fprintf(stderr,"Can`t multiply Fp-matrix and Z-matrix\n");
    exit(3);
  }
  T_mat = mat_mul(L_mat, M_mat);
  P_mat = mat_mul(T_mat, R_mat);
  free_mat(T_mat);
  return(P_mat);
}

/*}}}  */
