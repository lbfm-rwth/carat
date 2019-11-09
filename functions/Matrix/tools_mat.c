#include "typedef.h"
#include "tools.h"
#include "matrix.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: tools_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
/*{{{}}}*/
/*{{{  iset_entry*/
/*
@---------------------------------------------------------------------
@  boolean iset_entry(mat, r, c, v);
@ matrix_TYP *mat;
@ int r,c;
@ int v;
@
@  Set an integer entry in the r-th row and c-th column of mat to v.
@
@-------------------------------------------------------------------------
 */

boolean iset_entry(matrix_TYP *mat, int r, int c, int v)
{
int w;

  /* 
   * Check rows and columns 
   */
  if((r >= mat->rows) || (c >= mat->cols)) {
    fprintf(stderr,"Error: Limits of coordinates violated \n");
    return(FALSE);
  }

  if( mat->prime != 0 ) {
    if((w = v % mat->prime) < 0) {
      w += mat->prime;
    }
    mat->array.SZ[r][c] = w;
  } else if(mat->flags.Integral) {
    mat->array.SZ[r][c] = v;
  } else {
    mat->array.SZ[r][c] = v*mat->kgv;
  }
  if(mat->flags.Symmetric && (r != c)) { /* Check matrix flags */
    mat->flags.Symmetric = FALSE;
    Check_mat(mat);
  }
  return(TRUE);
}

/*}}}  */
/*{{{  rset_entry*/

/*
@-------------------------------------------------------------------------
@
@ boolean rset_entry(mat, r, c, v);
@ matrix_TYP *mat;
@ int r, c;
@ rational v;
@
@ Set an rational entry in the r-th row and c-th column of mat to v
@
@-------------------------------------------------------------------------
 */
boolean rset_entry(matrix_TYP *mat, int r, int c, rational v)
{
int t;
int **Z;
int w;
int i,j;

  if(v.n == 1) { /* rational is actually an integer */
    return(iset_entry(mat,r,c,v.z));
  }

  /* Check rows and columns */
  if((r >= mat->rows) || (c >= mat->cols)) {
    fprintf(stderr,"Error: Limits of coordinates violated \n");
    return(FALSE);
  }

  if( mat->prime != 0 ) { /* somehow pathologic, but why not */
    w= (v.z%v.n)%mat->prime;
    mat->array.SZ[r][c] = w;
  } else {
    t= GGT(mat->kgv,v.n);
    if(t == v.n) {
      mat->array.SZ[r][c] = v.z*mat->kgv/t;
    } else {
      int t1= v.n/t;
      mat->kgv *= t1;
      Z = mat->array.SZ;
      for(i = 0; i < mat->rows; i++) {
        for(j = 0; j < mat->cols; j++) {
          Z[i][j] *= t1;                
        }
      }
    Z[r][c]= v.z*mat->kgv/t;
    }
  }

  if(mat->flags.Symmetric && (r != c)) { /* Check matrix flags */
    mat->flags.Symmetric = FALSE;
    Check_mat(mat);
  }
  return(TRUE);
}
/*}}}  */
/*{{{  iscal_mul*/
/*
@-------------------------------------------------------------------------
@
@ void iscal_mul(mat, v);
@ matrix_TYP *mat;
@ int v;
@
@ Multiply matrix with integer scalar
@
@-------------------------------------------------------------------------
 */
void 
iscal_mul (matrix_TYP *mat, int v)
{
int i,j;
int **Z;
int t;

  /* handle the trivial cases */
  if (v == 1) {
    return;
  }
  if ( null_mat( mat ) ) {
    return;
  }
  if (v == 0) {
    t = 0;
    memset2dim( (char **)mat->array.SZ, mat->rows, mat->cols, sizeof(int), (char *)&t );
    Check_mat(mat);
    return;
  }

  t= mat->kgv%v;
  if ( (mat->prime == 0) && (t == 0)) {
    mat->kgv /= v;
    if(mat->kgv == 1 || mat->kgv == -1 ) { /* take care of negative values */
      mat->flags.Integral = TRUE;
      if(mat->kgv == -1) {
        mat->kgv = 1;
        Z = mat->array.SZ;
        if(mat->flags.Diagonal) {
          for(i = 0; i < mat->rows; i++) {
            Z[i][i] = -Z[i][i];           
          }
        } else {
          for(i = 0; i < mat->rows; i++) {
            for(j = 0; j < mat->cols; j++) {
              Z[i][j] = -Z[i][j];           
            }
          }
        }
      }
    }
    return;
  }
  
  /* now the real work */
  Z = mat->array.SZ;
  if( mat->prime != 0 ) {
    if((v %= mat->prime) < 0) {
      v += mat->prime;
    }
    init_prime(mat->prime);
    if(mat->flags.Diagonal) {
      for(i = 0; i < mat->rows; i++) {
        Z[i][i] = P(Z[i][i],v);       
      }
    } else {
      for(i = 0; i < mat->rows; i++) {
        for(j = 0; j < mat->cols; j++) {
          Z[i][j] = P(Z[i][j],v);
        }
      }
    }
  } else {  /* P_TYP */
    t= GGT(mat->kgv,v);
    if( t != 1 ) {
      mat->kgv /= t;
      v /= t;
    }
    if(mat->flags.Diagonal) {
      for(i = 0; i < mat->rows; i++) {
        Z[i][i] *= v;
      }
    } else {
      for(i = 0; i < mat->rows; i++) {
        for(j = 0; j < mat->cols; j++) {
          Z[i][j] *= v;                 
        }
      }
    }
  }
}

/*}}}  */
/*{{{  rscal_mul*/
/*                
@-------------------------------------------------------------------------
@
@ void rscal_mul(mat,v);
@ matrix_TYP *mat;
@ rational v;
@
@ Multiply matrix with rational scalar
@
@-------------------------------------------------------------------------
 */
void 
rscal_mul (matrix_TYP *mat, rational v)
{
int vz, t;
int pp;

  if( null_mat( mat ) ) {
    return;
  }
  if(v.z == 0) {
    t = 0;
    memset2dim( (char **)mat->array.SZ, mat->rows, mat->cols, sizeof(int), (char *)&t );
    Check_mat(mat);
    return;
  }
  t  = 0;
  Normal(&v);
  if(v.n == 1) {
    iscal_mul(mat, v.z);
    return;
  }
  
  if(mat->prime != 0) {
    init_prime(mat->prime);
    vz= v.z%mat->prime;
    t = v.n%mat->prime;
    pp = P(vz, t);
    iscal_mul(mat,pp);
  }
  if(mat->flags.Integral) {
    mat->kgv= v.n;
    mat->flags.Integral = FALSE;
    vz= v.z;
  } else {
    t= GGT(mat->kgv,v.z);
    mat->kgv= mat->kgv/t*v.n;
    vz= v.z/t;
  }
  
  if(v.z == 1) {
    Check_mat(mat);
  } else {
    iscal_mul(mat,vz);
  }
}
/*}}}  */
/*{{{  kill_row*/
/*
@-------------------------------------------------------------------------
@
@ boolean kill_row(mat, row);
@ matrix_TYP *mat;
@ int row;
@
@ Kill one row in the matrix
@ The argument 'row' must use the indexing of the matrix, ie.
@ the rows are numbered from 0 to n-1.
@
@-------------------------------------------------------------------------
 */
boolean kill_row(matrix_TYP *mat, int row)
{
int i;
int **Z;

  if(row >= mat->rows) {
    fprintf(stderr,"Error in kill_row: not so many rows!\n");
    return(FALSE);
  }
  
  Z = mat->array.SZ;
  free( Z[row] );
  for(i = row+1; i < mat->rows; i++) {
    Z[i-1] = Z[i];
  }                                  
  mat->rows--;
  mat->array.SZ = (int **)realloc(Z,mat->rows*sizeof(int  *));
  
  if( mat->array.N != NULL ) {
    Z = mat->array.N;
    free(Z[row]);
    for(i = row+1; i < mat->rows+1; i++) {
      Z[i-1] = Z[i];
    }
    mat->array.N = (int **)realloc(Z,mat->rows*sizeof(int *));
  }       
  mat->flags.Symmetric = FALSE;
  Check_mat(mat);
  return(TRUE);
}
/*}}}  */
/*{{{  kill_col*/
/*
@-------------------------------------------------------------------------
@
@ boolean kill_col(mat, col);
@ matrix_TYP *mat;
@ int col;
@
@ Kill one column in the matrix
@ The argument 'col' must use the indexing of the matrix, ie.
@ the cols are numbered from 0 to n-1.
@
@-------------------------------------------------------------------------
 */
boolean kill_col(matrix_TYP *mat, int col)
{
int i, j;
int **Z;

  if(col >= mat->cols) {
    fprintf(stderr,"Error in kill_col: not so many cols!\n");
    return(FALSE);
  }
  
  Z = mat->array.SZ;
  for(i = 0; i < mat->rows; i++) {
    for(j = col; j < mat->cols - 1; j++) {
      Z[i][j] = Z[i][j+1];
    }
    Z[i] = (int *)realloc(Z[i], (mat->cols-1)*sizeof(int));
  }
  
  if( mat->array.N != NULL) {
    Z = mat->array.N;
    for(i = 0; i < mat->rows; i++) {
      for(j = col; j < mat->cols - 1; j++) {
        Z[i][j] = Z[i][j+1];
      }
      Z[i] = (int *)realloc(Z[i], (mat->cols-1)*sizeof(int));
    }
  }
  
  mat->cols--;
  mat->flags.Symmetric = FALSE;
  Check_mat(mat);
  return(TRUE);
}

/*}}}  */
/*{{{  ins_row*/
/*
@-------------------------------------------------------------------------
@ boolean ins_row(mat, row);
@ matrix_TYP *mat;
@ int row;
@
@ Insert one row in the matrix
@ The argument 'row' must use the indexing of the matrix, ie.
@ the rows are numbered from 0 to n-1.
@
@-------------------------------------------------------------------------
 */
boolean ins_row(matrix_TYP *mat, int row)
{
int i;
int **Z;

  if(row > mat->rows) {
    fprintf(stderr,"Error in ins_row: not so many rows!\n");
    return(FALSE);
  }
  
  mat->rows++;
  Z=(int **)realloc(mat->array.SZ,mat->rows*sizeof(int *));
  mat->array.SZ = Z;
  for(i = mat->rows-1; i > row; i--) {
    Z[i] = Z[i-1];
  }
  Z[row] = (int  *)calloc(mat->cols , sizeof(int ));

  if( mat->array.N) {
    Z=(int **)realloc(mat->array.N,mat->rows*sizeof(int *));
    mat->array.N = Z;
    for(i = mat->rows-1; i > row; i--) {
      Z[i] = Z[i-1];
    }
    Z[row] = (int  *)calloc(mat->cols , sizeof(int ));
  }
  
  Check_mat(mat);
  return(TRUE);
}
/*}}}  */
/*{{{  ins_col*/
/*
@-------------------------------------------------------------------------
@ boolean ins_col(mat, col);
@ matrix_TYP *mat;
@ int col;
@
@ Insert one column in the matrix
@ The argument 'col' must use the indexing of the matrix, ie.
@ the cols are numbered from 0 to n-1.
@
@-------------------------------------------------------------------------
 */

boolean ins_col(matrix_TYP *mat, int col)
{
int i, j;
int **Z;

  if(col > mat->cols) {
    fprintf(stderr,"Error in kill_col: not so many cols!\n");
    return(FALSE);
  }
  
  mat->cols++;
  
  Z = mat->array.SZ;
  for(i = 0; i < mat->rows; i++) {
    Z[i] = (int *)realloc(Z[i], mat->cols*sizeof(int));
    for(j = mat->cols-1; j > col; j--) {
      Z[i][j] = Z[i][j-1];
    }
    Z[i][col] = 0;
  }
  if( mat->array.N != NULL) {
    Z = mat->array.N;
    for(i = 0; i < mat->rows; i++) {
      Z[i] = (int  *)realloc(Z[i], mat->cols*sizeof(int));
      for(j = mat->cols-1; j > col; j--) {
        Z[i][j] = Z[i][j-1];
      }
      Z[i][col] = 0;
    }
  }
  
  Check_mat(mat);
  return(TRUE);
}
/*}}}  */
/*{{{  imul_row*/
/*
@-------------------------------------------------------------------------
@ boolean imul_row(mat, row, v);
@ matrix_TYP *mat;
@ int row, v;
@
@ Multiply row with an integer
@
@-------------------------------------------------------------------------
 */
boolean imul_row(matrix_TYP *mat, int row, int v)
{
int j;
int **Z;

  if(row >= mat->rows) {
    fprintf(stderr,"Error in imul_row: not so many rows!\n");
    return(FALSE);
  }
  /* handle the trivial cases */
  if (v == 1) {
    return(TRUE);
  }
  if ( null_mat( mat ) ) {
    return(TRUE);
  }
  if (v == 0) {
    memset( mat->array.SZ[row], '\0', sizeof(int)*mat->cols );
    Check_mat(mat);
    return(TRUE);
  }
  Z = mat->array.SZ;
  if( mat->prime != 0) {
    if( (v %= mat->prime) < 0 ) {
      v += mat->prime;
    }
    if(v == 1) {
      return(TRUE);
    }
    init_prime(mat->prime);
    if(mat->flags.Diagonal) {
        Z[row][row] = P(Z[row][row],v);
    } else {
      for(j = 0; j < mat->cols; j++) {
        Z[row][j] = P(Z[row][j],v);   
      }
    }
  } else { /* flags & P_TYP */
    if(mat->flags.Diagonal) {
      Z[row][row] *= v;
    } else {
      for(j = 0; j < mat->cols; j++) {
        Z[row][j] *= v;
      }
    }
  }
  
  if(!(mat->flags.Diagonal)) {
    mat->flags.Symmetric = FALSE;
  }
  Check_mat(mat);
  return(TRUE);
}

/*}}}  */
/*{{{  rmul_row*/
/*
@-------------------------------------------------------------------------
@ boolean rmul_row(mat, row, v);
@ matrix_TYP *mat;
@ int row;
@ rational v;
@
@ Multiply row with a rational
@
@-------------------------------------------------------------------------
 */
boolean rmul_row(matrix_TYP *mat, int row, rational v)
{
int i,j;
int waste, t;
int **Z;

  if(row >= mat->rows) {
    fprintf(stderr,"Error in rmul_row: not so many rows!\n");
    return(FALSE);
  }
  if( mat->prime != 0 ) {
    t= (v.z%v.n)%mat->prime;
  /*  modeq(mod(&t,v.z,v.n),int_to_lang(mat->prime)); */
    return(imul_row(mat,row,t));
  }
  imul_row(mat,row,v.z);
  if(v.n == 1) {
    return(TRUE);
  }
  
  Z = mat->array.SZ;
  if(mat->flags.Diagonal) {
    if(Z[row][row] == 0) {
      return(TRUE);       
    }
    t= GGT(Z[row][row], v.n);
    if(t == v.n) {
      Z[row][row] /= v.n;
      return(TRUE);
    } else {
      waste= v.n/t;
      Z[row][row] /= t;
      mat->kgv *= waste;
      for(i = 0; i < row; i++) {
        Z[i][i] *= waste;       
      }
      for(i = row+1; i < mat->rows; i++) {
        Z[i][i] *= waste;                 
      }
    }
  } else {
    t= v.n;
    for(i = 0; i < mat->cols; i++) {
      if(Z[row][i] != 0) {
        t= GGT(t,Z[row][i]);
        if(t == 1) {
          i = mat->cols;
        }
      }
    }
    if(t != 1) {
      for(i = 0; i < mat->cols; i++) {
        if(Z[row][i] != 0) {
          Z[row][i] /= t;  
        }
      }
      waste= v.n/t;
    } else {
      waste= v.n;
    }
    if(waste != 1) {
      mat->kgv *= waste;
      for(j = 0; j < mat->cols; j++) {
        for(i = 0; i < row; i++) {
          if(Z[i][j] != 0) Z[i][j] *= waste;
        }
        for(i = row+1; i < mat->rows; i++) {
          if(Z[i][j] != 0) Z[i][j] *= waste;
        }
      }
    }
  }
  
  if(!(mat->flags.Diagonal)) {
    mat->flags.Symmetric = FALSE;
  }
  if(mat->kgv != 1) {
    mat->flags.Integral = FALSE;
  }
  Check_mat(mat);
  return(TRUE);
}
/*}}}  */
/*{{{  imulcol*/
/*
@-------------------------------------------------------------------------
@boolean imul_col(mat, col, v)
@ matrix_TYP *mat;
@ int col, v;
@
@ Multiply col with an integer
@
@-------------------------------------------------------------------------
 */
boolean imul_col(matrix_TYP *mat, int col, int v)
{
int j;
int **Z;

  if(col >= mat->cols) {
    fprintf(stderr,"Error in imul_col: not so many cols!\n");
    return(FALSE);
  }
  /* handle the trivial cases */
  if (v == 1) {
    return(TRUE);
  }
  if ( null_mat( mat ) ) {
    return(TRUE);
  }
  if (v == 0) {
    memset( mat->array.SZ[col], '\0', sizeof(int)*mat->rows );
    Check_mat(mat);
    return(TRUE);
  }

  Z = mat->array.SZ;
  if( mat->prime != 0 ) {
    if((v %= mat->prime) < 0) {
      v += mat->prime;
    }
    if(v == 1) {
      return(TRUE);
    }
    init_prime(mat->prime);
    if(mat->flags.Diagonal) {
        Z[col][col] = P(Z[col][col],v);
    } else {
      for(j = 0; j < mat->rows; j++) {
        Z[j][col] = P(Z[j][col],v);
      }
    }
  } else {
    if(mat->flags.Diagonal) {
      Z[col][col] *= v;
    } else {
      for(j = 0; j < mat->rows; j++) {
        Z[j][col] *= v;
      }
    }
  }
  if(!(mat->flags.Diagonal)) {
    mat->flags.Symmetric = FALSE;
  }
  Check_mat(mat);
  return(TRUE);
}

/*}}}  */
/*{{{  rmul_col*/
/*
@-------------------------------------------------------------------------
@ boolean rmul_col(mat, col, v);
@ matrix_TYP *mat;
@ int col;
@ rational v;
@
@ Multiply column with a rational
@
@-------------------------------------------------------------------------
 */
boolean rmul_col(matrix_TYP *mat, int col, rational v)
{
int i,j;
int waste, t;
int **Z;

  if(col >= mat->cols) {
    fprintf(stderr,"Error in rmul_col: not so many cols!\n");
    return(FALSE);
  }
  if( mat->prime != 0 ) {
    t= (v.z%v.n)%mat->prime;
    return(imul_col(mat,col,t));
  }
  
  imul_col(mat,col,v.z);
  if(v.n == 1) {
    return(TRUE);
  }
  
  Z = mat->array.SZ;
  if(mat->flags.Diagonal) {
    if(Z[col][col] == 0) {
      return(TRUE);      
    }
    t= GGT(Z[col][col], v.n);
    if(t == v.n) {
      Z[col][col] /= v.n;
      return(TRUE);
    } else {
      waste= v.n/t;
      Z[col][col] /= t;
      mat->kgv *= waste;
      for(i = 0; i < col; i++) {
        Z[i][i] *= waste;      
      }
      for(i = col+1; i < mat->cols; i++) {
        Z[i][i] *= waste;                
      }
    }
  } else {
    t= v.n;
    for(i = 0; i < mat->rows; i++) {
      if(Z[i][col] != 0) {
        t= GGT(t,Z[i][col]);
        if(t == 1) { 
          i = mat->cols;
        }
      }
    }
    if(t != 1) {
      for(i = 0; i < mat->rows; i++) {
        if(Z[i][col] != 0) {
          Z[i][col] /= t;  
        }
      }
      waste= v.n/t;
    } else {
      waste= v.n;
    }
    if(waste != 1) {
      for(i = 0; i < mat->rows; i++) {
        if(Z[i][col] != 0) {
          Z[i][col] /= waste;
        }
      }
      mat->kgv *= waste;
      for(j = 0; j < mat->rows; j++) {
        for(i = 0; i < col; i++) {
          Z[j][i] *= waste;      
        }
        for(i = col+1; i < mat->cols; i++) {
          Z[j][i] *= waste;                
        }
      }
    }
  }
  
  if(!(mat->flags.Diagonal)) {
    mat->flags.Symmetric = FALSE;
  }
  if(mat->kgv != 1) {
    mat->flags.Integral = FALSE;
  }
  Check_mat(mat);
  return(TRUE);
}
/*}}}  */

/*{{{  iadd_row*/
/*               
@-------------------------------------------------------------------------
@  boolean iadd_row(mat,t_row, d_row, v)
@ matrix_TYP *mat;
@ int t_row, d_row, v;
@
@  Add v-times t_row to d_row of mat, where v is an integer
@
@-------------------------------------------------------------------------
 */
boolean iadd_row(matrix_TYP *mat, int t_row, int d_row, int v)
{
int **Z, i,j,t;

  if((t_row >= mat->rows) || (d_row >= mat->rows)) {
    fprintf(stderr,"Error in iadd_row: not so many rows!\n");
    return(FALSE);
  }
  
  if( mat->prime ) {
    if((v %= mat->prime) < 0) v += mat->prime;
  }
  if(v == 0) {
    return(TRUE);
  }
  if((mat->prime == 0) && (mat->kgv == 0)) {
    return(TRUE);                          
  }
  
  if(d_row == t_row) {
    fprintf(stderr,"Warning: Same t_row and d_row in iadd_row\n");
    return(imul_row(mat,t_row,v+1));
  }
  
  Z = mat->array.SZ;
  if( mat->prime != 0 ) {
    init_prime(mat->prime);
    if(mat->flags.Diagonal) {
      Z[d_row][t_row] = P(v,Z[t_row][t_row]);
      if(Z[d_row][t_row]) {
        mat->flags.Symmetric =
        mat->flags.Diagonal = FALSE;
      }
    } else {
      for(i = 0; i < mat->cols; i++) {
        Z[d_row][i] = S(Z[d_row][i],P(v,Z[t_row][i]));
      }
    }
    Check_mat(mat);
    return(TRUE);
  }
  if(mat->flags.Diagonal) {
    if(Z[t_row][t_row] ) {
      i = Z[t_row][t_row] * v;
      Z[d_row][t_row] = i;
      mat->flags.Symmetric =
      mat->flags.Diagonal  = FALSE;
    }
  } else {  /* diagonal */
    j = find_max_entry( mat );
    for(i = 0; i < mat->cols; i++) {
      Z[d_row][i] += v * Z[t_row][i];
      t = abs(Z[d_row][i]);
      if(t > j) {
        j = t;
      }
    }
    mat->flags.Symmetric = FALSE;
  }
  Check_mat(mat);
  return(TRUE);
}
/*}}}  */
/*{{{  radd_row*/
/*
@-------------------------------------------------------------------------
@ boolean radd_row(mat,t_row, d_row, v)
@ matrix_TYP *mat;
@ int t_row, d_row;
@ rational v;
@
@ Add v-times t_row to d_row of mat, where v is a rational
@
@-------------------------------------------------------------------------
 */
boolean radd_row(matrix_TYP *mat, int t_row, int d_row, rational v)
{
int **Z, i, j;
int waste,t;

  if((t_row >= mat->rows) || (d_row >= mat->rows)) {
    fprintf(stderr,"Error in radd_row: not so many rows!\n");
    return(FALSE);
  }
  t = 0;
  
  if( mat->prime != 0) {
    waste= (v.z%v.n)%mat->prime;
    return(iadd_row(mat,t_row, d_row, waste));
  }
  if(v.z == 0) {
    return(TRUE);
  }
  if( null_mat( mat ) ) {
    return(TRUE);            
  }
  if(v.n == 1) {
    return(iadd_row(mat,t_row, d_row, v.z));
  }
  
  if(d_row == t_row) {
    fprintf(stderr,"Warning: Same t_row and d_row in ladd_row\n");
    v.z += v.n;
    i = rmul_row(mat,t_row,v);
    v.z -= v.n;
    return(i);
  }
  
  Z = mat->array.SZ;
  if(mat->flags.Diagonal) {
    if(Z[t_row][t_row] != 0) {
      mat->flags.Symmetric =
      mat->flags.Diagonal = FALSE;
      Z[d_row][t_row]= v.z*Z[t_row][t_row];
      waste= GGT(Z[d_row][t_row],v.n);
      if(v.n == waste) {
        Z[d_row][t_row] /= v.n;
      } else {
        t= v.n/waste;
        Z[d_row][t_row] /= waste;
        mat->kgv *= t;
        for(i = 0; i < mat->rows; i++) {
          Z[i][i] *= t;
        }
        mat->flags.Integral = FALSE;
      }
    }
  } else {
    t= v.n;
    for(i = 0; i < mat->cols; i++) {
      Z[d_row][i]=Z[d_row][i]*v.n+Z[t_row][i]*v.z;
      if((t != 1) && (Z[d_row][i] != 0)) {
        t= GGT(Z[d_row][i],t);
      }
    }
    if(t == v.n) {/* the lucky case: row-entries divisible by v.n! */
      for(i = 0; i < mat->cols; i++) {
        Z[d_row][i] /= v.n;          
      }
    } else if( t == 1) { /* the opposite case: no common factor */
      mat->kgv *= v.n;
      for(j = 0; j < mat->cols; j++) {
        for (i = 0; i < d_row; i++) {
          Z[i][j] *= v.n;           
        }
        for (i = d_row+1; i < mat->rows; i++) {
          Z[i][j] *= v.n;                     
        }
      }
      mat->flags.Integral = FALSE;
    } else { /* something between the last two cases: most work */
      waste = v.n / t;
      mat->kgv *=  waste;
      for(j = 0; j < mat->cols; j++) {
        for (i = 0; i < d_row; i++) {
          Z[i][j] *= waste;         
        }
        for (i = d_row+1; i < mat->rows; i++) {
          Z[i][j] *= waste;                   
        }
        Z[d_row][j] /= t;
      }
      mat->flags.Integral = FALSE;
    }
  }
  
  Check_mat(mat);
  return(TRUE);
}
/*}}}  */


/*{{{  iadd_col*/
/*
@-------------------------------------------------------------------------
@ boolean iadd_col(mat,t_col, d_col, v)
@ matrix_TYP *mat;
@ int t_col, d_col, v;
@
@ Add v-times t_col to d_col of mat, where v is an integer
@
@-------------------------------------------------------------------------
 */
boolean iadd_col(matrix_TYP *mat, int t_col, int d_col, int v)
{
int **Z, i,j,t;

  if((t_col >= mat->cols) || (d_col >= mat->cols)) {
    fprintf(stderr,"Error in iadd_row: not so many rows!\n");
    return(FALSE);
  }
  if( mat->prime != 0 ) {
    if((v %= mat->prime) < 0) v += mat->prime;
  }
  if(v == 0) {
    return(TRUE);
  }

  if( null_mat( mat ) ) {
    return(TRUE);                               
  }
  
  if(d_col == t_col) {
    fprintf(stderr,"Warning: Same t_col and d_col in iadd_col\n");
    return(imul_col(mat,t_col,v+1));
  }
  
  Z = mat->array.SZ;
  if( mat->prime != 0 ) {
    init_prime(mat->prime);
    if(mat->flags.Diagonal) {
      Z[t_col][d_col] = P(v,Z[t_col][t_col]);
      if(Z[t_col][d_col]) {
        mat->flags.Symmetric =
        mat->flags.Diagonal = FALSE;
      }
    } else {
      for(i = 0; i < mat->rows; i++) {
        Z[i][d_col] = S(Z[i][d_col],P(v,Z[i][t_col]));
      }
    }
    Check_mat(mat);
    return(TRUE);
  }
  if(mat->flags.Diagonal) {
    if(Z[t_col][t_col] ) {
      i = Z[t_col][t_col] * v;
      Z[t_col][d_col] = i;
      mat->flags.Symmetric =
      mat->flags.Diagonal = FALSE;
    }
  } else {
    j = find_max_entry( mat );
    for(i = 0; i < mat->rows; i++) {
      Z[i][d_col] += v * Z[i][t_col];
      t = abs(Z[i][d_col]);
      if(t > j) {
        j = t;
      }
    }
    mat->flags.Symmetric = FALSE;
  }
  Check_mat(mat);
  return(TRUE);
}
/*}}}  */
/*{{{  radd_col*/
/*
@-------------------------------------------------------------------------
@ boolean radd_col(mat,t_col, d_col, v)
@ matrix_TYP *mat;
@ int t_col, d_col;
@ rational v;
@
@  Add v-times t_col to d_col of mat, where v is a rational
@
@-------------------------------------------------------------------------
 */
boolean radd_col(matrix_TYP *mat, int t_col, int d_col, rational v)
{
int **Z, i,j;
int waste,t;

  if((t_col >= mat->cols) || (d_col >= mat->cols)) {
    fprintf(stderr,"Error in radd_col: not so many cols!\n");
    return(FALSE);
  }
  t = 0;
  
  if( mat->prime != 0) {
    return(iadd_col(mat,t_col, d_col, (v.z%v.n)%mat->prime));
  }
  if(v.z == 0) {
    return(TRUE);
  }
  if( null_mat( mat ) ) {
    return(TRUE);             
  }
  if(v.n == 1) {
    return(iadd_col(mat,t_col, d_col, v.z));
  }
  
  if(d_col == t_col) {
    fprintf(stderr,"Warning: Same t_col and d_col in radd_col\n");
    v.z += v.n;
    i = rmul_row(mat,t_col,v);
    v.z -= v.n;
    return(i);
  }
  
  Z = mat->array.SZ;
  if(mat->flags.Diagonal) {
    if(Z[t_col][t_col] != 0) {
      mat->flags.Symmetric =
      mat->flags.Diagonal = FALSE;
      Z[t_col][d_col] = v.z * Z[t_col][t_col];
      waste = GGT(Z[t_col][d_col],v.n);
      if(v.n == waste) {
        Z[t_col][d_col]/= v.n;
      } else {
        t = v.n / waste;
        Z[t_col][d_col] /= waste;
        mat->kgv *= t;
        for(i = 0; i < mat->rows; i++) {
          Z[i][i] *= t;
        }
        mat->flags.Integral = FALSE;
      }
    }
  } else {
    t = v.n;
    for(i = 0; i < mat->rows; i++) {
      Z[i][d_col] = Z[i][d_col] * v.n + Z[i][t_col] * v.z;
      if((t != 1) && (Z[i][d_col] != 0)) {
        t= GGT(Z[i][d_col],t);
      }
    }
    if( t == v.n ) {/* the lucky case: row-entries divisible by v.n! */
      for(i = 0; i < mat->rows; i++) {
        Z[i][d_col] /= v.n;          
      }
    } else if( t == 1) { /* the opposite case: no common factor */
      mat->kgv *= v.n;
      for(j = 0; j < mat->rows; j++) {
        for (i = 0; i < d_col; i++) {
          Z[j][i] *= v.n;           
        }
        for (i = d_col+1; i < mat->cols; i++) {
          Z[j][i] *= v.n;                     
        }
      }
      mat->flags.Integral = FALSE;
    } else { /* something between the last two cases: most work */
      waste= v.n / t;
      mat->kgv *= waste;
      for(j = 0; j < mat->rows; j++) {
        for (i = 0; i < d_col; i++) {
          Z[j][i] *= waste;         
        }
        for (i = d_col+1; i < mat->cols; i++) {
          Z[j][i] *= waste;                    
        }
        Z[j][d_col] /= t;
      }
      mat->flags.Integral = FALSE;
    }
  }
  Check_mat(mat);
  return(TRUE);
}
/*}}}  */
/*{{{  normal_rows*/
/* 
@-------------------------------------------------------------------------
@ void normal_rows(mat);
@ matrix_TYP *mat;
@
@ Tool for integral Gauss-elimination: the entries of the
@ rows of mat are divided by their row-gcd.
@
@-------------------------------------------------------------------------
 */
void 
normal_rows (matrix_TYP *mat)
{
int i,j, g;
int **Z, *S_I;

  if( mat->flags.Diagonal || (mat->prime != 0) || null_mat( mat) ) {
    return;                                                          
  }
  Z = mat->array.SZ;
  for(i = 0; i < mat->rows; i++) {
    S_I = Z[i];
    g = 0;
    for(j = 0; j < mat->cols; j++) {
      if(S_I[j]) g = GGT(g,S_I[j]);
      if(g == 1) j = mat->cols;
    }
    if(g != 1) {
      for(j = 0; j < mat->cols; j++) {
        if(S_I[j]) {
          S_I[j] /= g;
        }
      }
    }
  }    
}
/*}}}  */
/*{{{  normal_cols*/
/* 
@-------------------------------------------------------------------------
@ void normal_cols(mat)
@
@ Analogue to normal_rows() function for cols
@
@-------------------------------------------------------------------------
 */
void 
normal_cols (matrix_TYP *mat)
{
int i,j, g;
int **Z;

  if( mat->flags.Diagonal || (mat->prime != 0) || null_mat( mat ) ) {
    return;                                                                 
  }

  Z = mat->array.SZ;
  for(i = 0; i < mat->cols; i++) {
    g = 0;
    for(j = 0; j < mat->rows; j++) {
      if(Z[j][i]) {
        g = GGT(g,Z[j][i]);
      }
      if(g == 1) {
        j = mat->rows;
      }
    }
    if(g != 1) {
      for(j = 0; j < mat->rows; j++) {
        if(Z[j][i]) { 
          Z[j][i] /= g;
        }
      }
    }
  }
}
/*}}}  */
/*{{{  mat_to_line*/
/*
@-------------------------------------------------------------------------
@ matrix_TYP *mat_to_line(gen, num);
@ matrix_TYP **gen;
@ int num;
@ 
@ matrix_TYP **gen;
@ int num;
@
@ converts each matrix of gen into a row-vector, forgets the kgv.
@ the i-th matrix is converted to the i-th row of the result.
@ Tool for determine a basis of a vector space of matrices
@
@-------------------------------------------------------------------------
 */
matrix_TYP *
mat_to_line (matrix_TYP **gen, int num)
{
matrix_TYP *mat;
int i, j, k, row, col;
int **Z, **ZI;

  row = gen[0]->rows;
  col = gen[0]->cols;
  for(i = 0; i < num; i++) {
    if((row != gen[i]->rows) || (col != gen[i]->cols)) {
      fprintf(stderr,"Fehler: verschiedene Dimensionen in mat_to_line\n");
      exit(3);
    }
  }
  mat = init_mat(num, row * col,"");
  Z = mat->array.SZ;
  for(k = 0; k < num; k++) {
    ZI = gen[k]->array.SZ;
    for(j = 0; j < row; j++) {
      memcpy(Z[k]+j*col,ZI[j],col*sizeof(int));
    }
  }
  return(mat);
}

/*}}}  */
/*{{{  line_to_mat*/
/*
@-------------------------------------------------------------------------
@ matrix_TYP **line_to_mat(mat,row,col)
@ matrix_TYP *mat;
@ int row, col;
@
@ Inverse function to mat_to_line
@
@-------------------------------------------------------------------------
 */
matrix_TYP **
line_to_mat (matrix_TYP *mat, int row, int col)
{
matrix_TYP **gen;
int i,k;
int **Z, **ZI;

  if(mat->cols != row* col) {
    fprintf(stderr,"Fehler in line_to_mat: Falsche Dimensionierung\n");
    exit(3);
  }
  gen = (matrix_TYP **)malloc(mat->rows * sizeof(matrix_TYP *));
  
  for(k = 0; k < mat->rows; k++) {
    Z = mat->array.SZ;
    gen[k] = init_mat(row, col, "");
    ZI = gen[k]->array.SZ;
    for(i = 0; i < mat->cols; i++) {
      ZI[i / col][i % col] = Z[k][i];
    }
    Check_mat(gen[k]);
  }
  return(gen);
}

/*}}}  */
/*{{{  normal_mat*/
/*
@-------------------------------------------------------------------------
@ int normal_mat(mat);
@ matrix_TYP *mat;
@
@ Calculate ggt of mat-entries and divide matrix by it
@ Return value is the ggt
@
@-------------------------------------------------------------------------
 */
int 
normal_mat (matrix_TYP *mat)
{
int i,j, g;
int  **Z, *S_I;

  if( null_mat(mat) || (mat->prime != 0)) {
    return(1);
  }

  Z = mat->array.SZ;
  g = 0;
  if(mat->flags.Diagonal) {
    for(i = 0; i < mat->rows; i++) {
      if(Z[i][i]) g = GGT(g,Z[i][i]);
      if(g == 1) {
        i = mat->rows;
      }
    }
  } else {
    g = 0;
    for(i = 0; i < mat->rows; i++) {
      S_I = Z[i];
      for(j = 0; j < mat->cols; j++) {
        if(S_I[j]) {
          g = GGT(g,S_I[j]);
        }
        if(g == 1) {
          j = mat->cols;
          i = mat->rows;
        }
      }
    }
  }
  if(g != 1) {
    for(i = 0; i < mat->rows; i++) {
      S_I = Z[i];
      for(j = 0; j < mat->cols; j++) {
        if(S_I[j]) S_I[j] /= g;      
      }
    }
  }
  return(g);
}       
/*}}}  */

