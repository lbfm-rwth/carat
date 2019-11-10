#include "typedef.h"
#include "matrix.h"
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: construct_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
/*
 |
 | matrices/construct.c -- Funktionen, die Matrizen erzeugen
 |
 |
 | exportiert die Funktionen
 |
 | init_mat
 | copy_mat
 | free_mat
 | Check_mat
 |
 | importiert die Funktionen
 | GGT, calloc2dim, malloc2dim, memcpy2dim, free2dim
 | importiert die Variable
 | act_prime
 |
*/

/*{{{}}}*/
/*{{{  init_mat*/
/*{{{  Documentation*/
/*
 @
 @ void init_mat( cols, rows, Optionen );
 @ generaates a matrix with the options described below
 | Erzeugt eine Matrix gemaess der uebergebenen Parameter, s.u.
 @
 @ Args:
 @ int cols, int rows -- number of rows and columns
 @ char *Optionen     -- String consistin of numbers and letters 'c', 'e', 'd',
 @                       's', 'r', 'p'
 | Die Bedeutung der einzelnen Optionen:
 @ Explanation of the options:
 @
 @                Options
 @  ==========================================================================
 @  's'     : symmetric matrix                    |  mat->flags.Symmetric=1
 @  --------------------------------------------------------------------------
 @  'd'     : diagonal matrix, implies 's'        |  mat->flags.Diagonal=1
 @                                                |  mat->flags.Symmetric=1
 @  --------------------------------------------------------------------------
 @  'c', 'e': scalar matrix, implies 'd' und 's'  |  mat->flags.Scalar=1
 @                                                |  mat->flags.Diagonal=1
 @                                                |  mat->flags.Symmetric=1
 @  --------------------------------------------------------------------------
 @     numerical letter from   0 .. 9.            |  mat->flags.Scalar=1
 @   The same as 'c' and 'e', but the matrix      |  mat->flags.Diagonal=1
 @   is initialized with the number               |  mat->flags.Symmetric=1
 @  --------------------------------------------------------------------------
 @  'p'     : Matrix over  GF(p) | mat->flags.Integral= 1
 @                               | mat->prime= act_prime (glob. Variable)
 @                               | mat->array.N= NULL
 @  --------------------------------------------------------------------------
 @  'r': rational  matrix,      | mat->flags.Integral= 0
 @                              | mat->array.N= malloc( was auch immer );
 @
 |  falls 'r' und 'p' gleichzeitig gesetzt werden, so hat 'p' Vorrang.
 |  mat->flags.Integral ist immer gesetzt, es sei denn, 'r' waere
 |  angegeben worden oder mat->kgv != 1
 @  if 'r' and 'p' are used at the same time, 'p' has priority
 @
 |  mat->kgv wird bei einer Bruchmatrix als Hauptnenner aller Eintraege
 |  benutzt. Bei einer Matrix ueber Z ist mat->kgv= 1. Falls mat->kgv != 1,
 |  so ist mat->array.N = NULL.
 @  mat->kgv is used by a rational matrix a least common divisor
 @  over the denomonatop of all entries. An integral matrix has mat->kgv = 1.
 @  If mat->kgv != 1 then mat->array.N = NULL.
 @
 |  Nullmatrizen werden nicht besonders gekennzeichnet. Das verschwendet
 |  zwar Rechenzeit, ist aber Programmtechnisch weniger aufwendig.
 |  Insgesamt sollte bei 6x6 Matrizen die Performance nicht allzusehr
 |  leiden.
 @  Zero matrices have no extra labeling
 @
 |  Wichtig: falls im Optionen-String keine Zahlen vorkommen, so erzeugt
 |  init_mat() eine Nullmatrix, alle Eintraege sind also mit '0' initialisiert
 @  if no numerical numbers are conatained in the string, all entries 
 @  of the matrix are initialized with 0.
 @
*/
/*}}}  */
/*{{{  Source-Code*/
matrix_TYP *
init_mat (int rows, int cols, const char *option)
{
int i;
matrix_TYP *mat;
int val;
const char *temp = 0;
const char *tempN;

  mat = NULL;
  mat = (matrix_TYP *)malloc(sizeof(matrix_TYP));

  /*{{{  parse the options*/
  if (strpbrk(option,"pP") != NULL)
  {
    mat->prime = act_prime; /* ist auf -1 gesetzt vor dem ersten Aufruf von */
                            /* init_prime().                                */
    mat->flags.Integral = 1;
  }
  else 
  { 
    mat->prime= 0;
    mat->flags.Integral  = strpbrk(option, "rR") == NULL ? TRUE: FALSE;
  }
  /* changed tilman 29/07/97 to cope with negative integers from
  if (    (strpbrk(option, "cC") != NULL) || (strpbrk(option, "eE") != NULL)
       || ( (temp = strpbrk(option, "0123456789")) != NULL )
     ) to : */
  if (    (strpbrk(option, "cC") != NULL) || (strpbrk(option, "eE") != NULL)
       || ( (temp = strpbrk(option, "-0123456789")) != NULL )
     )
  {
    mat->flags.Scalar=
    mat->flags.Diagonal=
    mat->flags.Symmetric= TRUE;
  }
  else
  {
    mat->flags.Scalar= FALSE;
    if ( strpbrk(option, "dD") != NULL )
    {
      mat->flags.Diagonal=
      mat->flags.Symmetric= TRUE;
    } 
    else
    {
      mat->flags.Diagonal= FALSE;
      mat->flags.Symmetric = strpbrk(option, "sS") != NULL ? TRUE: FALSE;
    }
  }
  if ( ( tempN = strpbrk(option, "/")) != NULL ) {
    if ( ( tempN = strpbrk(tempN, "0123456789")) != NULL ) {
      mat->flags.Integral = FALSE;
    }
  }
  /*}}}  */

  mat->kgv= 1;
  mat->cols = cols;
  mat->rows = rows;
  
  /*{{{  alloc SZ (is always done)*/
  mat->array.SZ= (int **) calloc2dim( rows, cols, sizeof(int) );
  if(temp != NULL) 
  {
    sscanf(temp,"%d", &val);
    for(i = 0; i < rows; i++) mat->array.SZ[i][i] = val;
  }
  /*}}}  */
  if ( !mat->flags.Integral )
  { 
    if (tempN != NULL) {
      sscanf(tempN,"%d", &val);
    } else {
      val = 1;
    }
    mat->array.N= (int **) calloc2dim( rows, cols, sizeof(int) );
    memset2dim( (char **)mat->array.N,rows,cols,sizeof(int), (char *)&val);
  }
  else 
    mat->array.N = NULL;
  
  return(mat);
}

/*}}}  */
/*}}}  */
/*{{{  copy_mat*/
/*
 @-----------------------------------------------------------------
 @ matrix_TYP *copy_mat( matrix_TYP *old );
 @
 | kopiert die Matrix "old", Rueckgabewert ist die Kopie
 @ copies the matrix 'old' and returns the copy.
 @-----------------------------------------------------------------
*/

matrix_TYP *
copy_mat (matrix_TYP *old)
{ 
matrix_TYP *newMat= NULL;

  if ( old )
  { 
    newMat= (matrix_TYP *)malloc( sizeof(matrix_TYP));
    if ( newMat )
    { 
      memcpy( (char *)newMat, (char *)old, sizeof(matrix_TYP) );
      if ( newMat->array.SZ )
      { 
        newMat->array.SZ= (int **)malloc2dim( old->rows, old->cols, sizeof(int) );
        memcpy2dim( (char **)newMat->array.SZ, (const char **)old->array.SZ,
                    old->rows, old->cols, sizeof(int) );
      }
      if ( newMat->array.N )
      { 
        newMat->array.N= (int **)malloc2dim( old->rows, old->cols, sizeof(int) );
        memcpy2dim( (char **)newMat->array.N,(const char **)old->array.N,
                    old->rows, old->cols, sizeof(int) );
      }
    }
  }
  return newMat;
}

/*}}}  */
/*{{{  free_mat*/
/*
 @--------------------------------------------------------------
 @ void free_mat ( mat )
 @ matrix_TYP *mat;
 | -- gibt den Speicher frei, den mat benutzt
 @
 @ clear the storage allocated for mat
 @--------------------------------------------------------------
*/

void 
free_mat (matrix_TYP *mat)
{  

  if ( mat->array.N )
    free2dim( (char **)mat->array.N, mat->rows);
  if ( mat->array.SZ )
    free2dim( (char **)mat->array.SZ, mat->rows);
  free( (int *)mat );

}

/*}}}  */
/*{{{  Check_mat*/
/*
 @--------------------------------------------------------------
 @ void Check_mat( mat )
 @ matrix_TYP *mat;
 @
 | Ueberprueft die mat->flags und korrigiert diese ggf. Kuerzt
 | Bruchmatrizen (sowohl matrix->array.N als auch mat->kgv )
 | Normiert die Darstellung modulo mat->prime, so dass 0 <= Eintrag < prime
 @ Checks mat->flags and corrects it if necessary.
 @ Reduces rational matrices (mat->array.N and also mat->kgv)
 @ Normalizes modulo mat->prime, such that 0<= entry < prime.
 @--------------------------------------------------------------
*/

void 
Check_mat (matrix_TYP *mat)
{
int i,j;
int g;
int **Z;

  if (mat->kgv == 0) mat->kgv = 1;
  if ( mat->array.N )
  {     
    if ( mat->kgv != 1 ) { 
      fprintf(stderr,"Check_mat: Error: denominator matrix with lcd != 1 detected.\n");
      exit(3);
    }      
    /*
     *  Beim Kuerzen wird Integral richtig gesetzt.
     */
    mat->flags.Integral= TRUE;
    /*{{{  kuerzen*/
    for ( i=0; i < mat->rows; i++ ) {
      for ( j=0; j < mat->cols; j++ ) {
        if ( mat->array.N[i][j] == 0 ) {
          fprintf (stderr,"Check_mat: Error: divide by zero\n");
          exit (3);
        }
        if ( mat->array.SZ[i][j] == 0 ) {
          mat->array.N[i][j]= 1;        
        } else {
          g = GGT (mat->array.SZ[i][j], mat->array.N[i][j]);
          mat->array.SZ[i][j] /= g;
          mat->array.N[i][j] /= g;
          if ( mat->array.N[i][j] < 0) {
            mat->array.N[i][j] = -mat->array.N[i][j];
            mat->array.SZ[i][j] = -mat->array.SZ[i][j];
          }
        }
        if ( mat->array.N[i][j] != 1 ) {
          mat->flags.Integral = FALSE;
        }
      }
    }
    /*}}}  */
    /*{{{  mat->flags setzen*/
    if ( mat->cols != mat->rows ) {
      mat->flags.Symmetric= mat->flags.Diagonal= mat->flags.Scalar= FALSE;
    } else {
      mat->flags.Diagonal= mat->flags.Symmetric= TRUE;
      for ( i=0; i < mat->rows-1 && mat->flags.Symmetric;i++ ) {
        for ( j=i+1;j < mat->cols && mat->flags.Symmetric; j++ ) {
          if ( mat->array.SZ[i][j] != mat->array.SZ[j][i] || 
               mat->array.N [i][j] != mat->array.N [j][i]    ) {
            mat->flags.Symmetric= mat->flags.Diagonal= FALSE;
          } else if ( mat->flags.Diagonal ) {
            mat->flags.Diagonal= mat->array.SZ[i][j] == 0;
          }
        }
      }
      mat->flags.Scalar = mat->flags.Diagonal;
      for ( i=0;i < mat->rows-1 && mat->flags.Scalar; i++) {
        if (    mat->array.SZ[i][i] != mat->array.SZ[i+1][i+1]
             || mat->array.N [i][i] != mat->array.N [i+1][i+1]   ) {
          mat->flags.Scalar= FALSE;                 
        }
      }
    } 
    /*
     *   release array.N if all denominators are equal to 1
     */
    if ( mat->flags.Integral ) {
      free2dim( (char **)mat->array.N, mat->rows );
      mat->array.N = NULL;
    }
    /*}}}  */
  } else { /* hoechstens noch Hauptnenner */
    Z = mat->array.SZ;
    if ( mat->prime != 0 ) { /* Matrix ueber GF(p) */
      /*{{{  normalisieren 0 <= zahl < p*/
      mat->flags.Integral= 1;
      if ( mat->prime != -1 ) /* passiert, wenn init_mat( ... , "p" ) vor */
      {                       /* init_prime() aufgerufen wird             */
        for (i = 0; i < mat->rows; i++)
          for (j = 0; j < mat->cols; j++)
            if ( (Z[i][j] %= mat->prime) < 0)
              Z[i][j] += mat->prime;
      }   
      /*}}}  */
    } 
    /*{{{  flags setzen*/
    if ( mat->cols != mat->rows )
      mat->flags.Symmetric= mat->flags.Diagonal= mat->flags.Scalar= FALSE;
    else
    {
      mat->flags.Diagonal= mat->flags.Symmetric= TRUE;
      for ( i=0; i < mat->rows-1 && mat->flags.Symmetric;i++ )
        for ( j=i+1;j < mat->cols && mat->flags.Symmetric; j++ )
        {
          if ( mat->array.SZ[i][j] != mat->array.SZ[j][i] )
          { 
            mat->flags.Symmetric= mat->flags.Diagonal= FALSE;
          } 
          else if ( mat->flags.Diagonal )
            mat->flags.Diagonal= mat->array.SZ[i][j] == 0;
        }
      mat->flags.Scalar = mat->flags.Diagonal;
      for ( i=0;i < mat->rows-1 && mat->flags.Scalar; i++)
      {
        if ( mat->array.SZ[i][i] != mat->array.SZ[i+1][i+1] )
          mat->flags.Scalar= FALSE;
      }
    }
    /*}}}  */
    if ( mat->prime == 0 )
    {
      /*{{{  kgv > 0 machen*/
      if ( mat->kgv < 0 )
      {
        mat->kgv= - mat->kgv;
        Z = mat->array.SZ;
        if(mat->flags.Diagonal)
          for (i = 0; i < mat->rows; i++) Z[i][i] = -Z[i][i];
        else
          for (i = 0; i < mat->rows; i++)
            for (j = 0; j < mat->cols; j++)
              Z[i][j] = -Z[i][j];
      }
      /*}}}  */
      /*{{{  kgv gegen alle Matrixelemente kuerzen*/
      g= mat->kgv;
      if( mat->flags.Scalar == FALSE )
      {
        if ( mat->flags.Diagonal == FALSE )
        {
          if ( mat->flags.Symmetric == FALSE )
          {
            for (i=0; i < mat->rows && g != 1;i++ )
              for ( j=0;j < mat->cols && g != 1; j++ )
                g= GGT ( g, mat->array.SZ[i][j] );
          }
          else
          {
            for (i=0; i < mat->rows && g != 1;i++ )
              for ( j=0;j <=i  && g != 1; j++ )
                g= GGT ( g, mat->array.SZ[i][j] );
          }
        }
        else
         for ( i=0; i < mat->rows && g != 1;i++) g= GGT( g, mat->array.SZ[i][i] );
      }
      /* timan 11/05/99: changed from
      else
      to : (to handle this rare case of a matrix without entries) */
      else if (mat->rows > 0 && mat->cols > 0)
        g = GGT(g, mat->array.SZ[0][0]);
      mat->kgv /= g;
      if(g != 1)
      {
         if ( mat->flags.Diagonal )
           for ( i=0; i < mat->rows;i++) mat->array.SZ[i][i] /= g;
         else
           for (i=0; i < mat->rows;i++ )
             for ( j=0;j < mat->cols; j++ )
               mat->array.SZ[i][j] /= g;
       }
       mat->flags.Integral = (mat->kgv == 1);
    }
  }
}

/*}}}  */

