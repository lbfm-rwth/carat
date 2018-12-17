/*{{{}}}*/
/*{{{  include*/
#include "typedef.h"
#include "tools.h"
#include "matrix.h"
#include "getput.h"
#include <string.h>

/*}}}  */

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: get_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
/*-----------------------------------------------------------*\
| Read matrix from infile                                     |
\*-----------------------------------------------------------*/

/*{{{  static alloc_N_if_neccessary*/
static void
alloc_N_if_neccessary( mat, flags )
matrix_TYP *mat;
flag_TYP *flags;
{

  if ( mat->array.N == NULL ) {
    if ( mat->prime != 0 ) {
      fprintf(stderr,"get_mat: Error: matrix over GF(%d) tries to be rational\n",
                     mat->prime);
      exit( 3 );
    } else if ( mat->kgv != 1 ) {
      fprintf(stderr,"get_mat: Error: matrix tries to use kgv %d and rational representation.\n",
                     mat->kgv );
      exit( 3 );
    } else {
      flags->Integral = FALSE;
      mat->array.N = (int **)malloc2dim( mat->rows, mat->cols,
                                         sizeof(int) );
      memset2dim( (char **) mat->array.N,
                mat->rows, mat->cols, sizeof(int), (char *)&Zero.n );
    }
  } 
}

/*}}}  */
/*{{{  fget_mat*/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *fget_mat (infile)
@ FILE *infile;
@ reads a matrix from infile
@---------------------------------------------------------------------------
@
\**************************************************************************/

matrix_TYP *fget_mat (infile)
FILE *infile;

{  
int rM, cM;
int **M;
matrix_TYP *mat;
flag_TYP flags;
char  string[256], *str ;
 int i, j, c;         




  flags.Integral  =
  flags.Symmetric =
  flags.Diagonal  = 
  flags.Scalar    = FALSE;
  /*------------------------------------------------------------*\
  | Read and scan header line                    |
  \*------------------------------------------------------------*/
  c=fscanf (infile, "%[ \t\n]", string);
  c=fscanf (infile, "%[^\n]",string);
  strtok (string, "%");
  sscanf (string, "%d", &rM);
  if ( (str = strpbrk (string, "xd")) != NULL ) {
    sscanf ( ++str, "%d", &cM);
    if ( str[-1] == 'x' ) {
      if ( cM == 0 ) {
        cM = rM;
        flags.Symmetric = TRUE;
      }
    } else {
      flags.Symmetric = flags.Diagonal = TRUE;
      flags.Scalar = (cM == 0);
      cM = rM;
    }
  }
  else {
    cM = rM;
  }
  /*
   * WARNING: You can read matrizes either of prime_Typ or with
   * kgv. If both are set in the headline, the second entry
   * will be ignored.
   */ 
  mat = init_mat(rM,cM,"");
  if ( (str = strchr (string, '/')) != NULL ) {
    if(str[1] == 'p') {
      str += 2;
      sscanf(str, "%d", &mat->prime);
      flags.Integral = TRUE;
    } else {
      sscanf(++str, "%d", &mat->kgv);
      flags.Integral = FALSE;
    }
  } else {
    flags.Integral = TRUE;  
  }    
  M = mat->array.SZ;
  for ( i = 0; i < rM; i++) {
    /*
     *  Read the matrix                       
     */
    if ( flags.Scalar ) {
      c=fscanf( infile, "%s", string );
      if ( strchr( string, '/') == NULL ) {
        sscanf(string, "%d", &M[0][0]);
        for ( i = 1; i < rM; i++ )M[i][i]= M[0][0];
      } else {                   
        alloc_N_if_neccessary( mat, &flags );
        sscanf(string, "%d/%d", &M[0][0],&mat->array.N[0][0]);
        for ( i = 1; i < rM; i++ ) {
          M[i][i]= M[0][0];
          mat->array.N[i][i] = mat->array.N[0][0];
        }
      }
    } else {
      if ( flags.Diagonal ) {
        for ( i = 0; i < cM; i++ ) {
          c=fscanf( infile, "%s", string );
          if ( strchr( string, '/') == NULL ) {
            sscanf(string, "%d", &M[i][i]);
          } else {                   
            alloc_N_if_neccessary( mat, &flags );
            sscanf(string, "%d/%d", &M[i][i],&mat->array.N[i][i]);
          }
        }
      } else {
        for ( i = 0; i < rM; i++  ) {
          for ( j = 0; j < ( flags.Symmetric ? i+1 : cM ); j++ ) {
            c=fscanf( infile, "%s", string );
            if ( strchr( string, '/') == NULL ) {
              sscanf(string, "%d", &M[i][j]);
            } else {                   
              alloc_N_if_neccessary( mat, &flags );
              sscanf(string, "%d/%d", &M[i][j],&mat->array.N[i][j]);
            }
          }
        }
        if ( flags.Symmetric ) {
          if ( mat->array.N != NULL ) {
            for ( i = 0; i < rM; i++) {
              for ( j = 0; j < i; j++) {
                M[j][i]= M[i][j];
                mat->array.N[j][i] = mat->array.N[i][j];
              }
            }
          } else {
            for ( i = 0; i < rM; i++) {
              for ( j = 0; j < i; j++) {
                M[j][i]= M[i][j];
              }
            }    
          }
        }
      }
    }
  }
  mat->flags.Integral  = flags.Integral;
  mat->flags.Symmetric = flags.Symmetric;
  mat->flags.Diagonal  = flags.Diagonal;
  mat->flags.Scalar    = flags.Scalar;
  Check_mat(mat);
  return ( mat );
}


/*}}}  */
/*{{{  get_mat*/

/**************************************************************************\
@---------------------------------------------------------------------------
@ matrix_TYP *get_mat (file_name)
@ char *file_name;
@ reads a matrix from the file with name 'file_name'
@---------------------------------------------------------------------------
@
\**************************************************************************/
matrix_TYP *get_mat (file_name)
char *file_name;

{  
matrix_TYP *mat;
FILE *infile;

/*------------------------------------------------------------*\
| Open input file                       |
\*------------------------------------------------------------*/
if ( file_name == NULL )
  infile = stdin;
else
  if ( (infile = fopen (file_name, "r")) == NULL ) {

    fprintf (stderr, "Could not open input-file %s\n", file_name);
    exit (4);
    }
/*------------------------------------------------------------*\
|  Call fget_mat                        |
\*------------------------------------------------------------*/
mat = fget_mat(infile);
  /*------------------------------------------------------------*\
  | Close input file                        |
  \*------------------------------------------------------------*/
if ( infile != stdin )
  fclose (infile);
return ( mat );
}

/*}}}  */
/*{{{  mget_mat*/
/*
@-------------------------------------------------------------------------
@
@  mat_array = mget_mat( file_name, anz );
@
@  differs from fmget_mat() only in specifying the file name instead of the
@  file-pointer.
@
@-------------------------------------------------------------------------
 */
matrix_TYP **mget_mat (file_name, anz)
char *file_name;
int *anz;
{  
matrix_TYP **mat;
FILE *infile;




  /*
   *    Open input file
   */
  if ( file_name == NULL ) {
    infile = stdin;        
  } else {
    if ( (infile = fopen (file_name, "r")) == NULL ) {
      perror("mget_mat: Error: Could not open input-file");
      exit (4);
    }
  }   
  mat = fmget_mat( infile, anz );
  /*
   * Close input file                        
   */
  if ( infile != stdin ) {
    fclose (infile);      
  }
  return ( mat );
}

/*}}}  */
/*{{{  fmget_mat*/
/*
@-------------------------------------------------------------------------
@
@ mat_array = fmget_mat( infile, anz );
@
@ matrix_TYP **mat_array: array of matrices read from infile
@
@ FILE *infile: the file the matrices are read from
@ int *anz:     the number of matrices that have been read from infile.
@               is set by the function.
@-------------------------------------------------------------------------
@
 */
matrix_TYP **fmget_mat (infile, anz)
FILE *infile;
int *anz;
{  
matrix_TYP **mat;
char string[512];
 int k, c;

  /*
   *   Open input file                       
   */
  c=fscanf (infile, "%[^\n]",string);
  if ( string[0] != '#' ) {
    *anz = 1;
    mat = (matrix_TYP **)malloc(sizeof(matrix_TYP *));
    rewind(infile);
    mat[0] = fget_mat(infile);
  } else {
    sscanf (string, "#%u", anz);
    mat = (matrix_TYP **)malloc(*anz*sizeof(matrix_TYP *));
    for ( k = 0; k < *anz; k++) {
      mat[k] = fget_mat(infile);
    }
  }
  return ( mat );
}
/*}}}  */
