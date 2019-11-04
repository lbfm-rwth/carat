#include "typedef.h"
#include "matrix.h" 
#include "tools.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: kgv2rat.
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
|
| kgv2rat.c
|
| exportiert die beiden Funktionen 
| result = kgv2rat( matrix );
| und
| result = rat2kgv( matrix );
|
| Die beiden Funktionen wandeln die Darstellung einer rationalen Matrix
| mittels des Hauptnenners mat->kgv und die Darstellung mittels der
| Nennermatrix mat->array.N ineinander um. 
| Das Ergebnis ist jeweils "maximal gekuerzt"
|
| Leider hat man in C keine Moeglichkeit, einen Integer-ueberlauf
| abzufangen, sonst koennte man in so einem Fall eine kgv-Matrix
| automatisch in eine Bruchmatrix umwandeln.
| 
|
\**************************************************************************/

/**************************************************************************\
@---------------------------------------------------------------------------
@ result = kgv2rat( matrix_TYP *mat );
@
@ changes an integral matrix or a matrix with mat->kgv != 1
@ to an rational matrix with allocated mat->array.N
@ the integer 'result' is 0, if this worked,
@ -1, if no storage for the matrix could be alloceted
@ -2, if mat->prime != 0
| wandelt eine Integer - oder kgv-Matrix in eine rational Matrix um
| "result" ist 0, falls alles glattging und nimmt den Wert
| -1 an, falls kein Speicher zur
| Erzeugung der Matrix vorhanden war.
| -2, falls mat->prime != 0
|
|
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
kgv2rat (matrix_TYP *mat)
{            
int i, j, result;

  if ( mat->prime ) /* this would be quite idiotic */
  { 
    result = -2;
  } 
  else if ( mat->array.N != NULL )
  {
/*
 * really nothing to do
 */
    mat->flags.Integral = FALSE;
    result = 0;
  } 
  else
  {      
    mat->array.N = (int **)malloc2dim( mat->rows, mat->cols, sizeof(int) );
    if ( mat->array.N == NULL )
    { 
      result = -1;

    } 
    else
    {
      mat->flags.Integral = FALSE;
      if ( mat->kgv == 0 )mat->kgv = 1;
      if ( mat->kgv == 1 )
      {
        for ( i=0; i < mat->rows; i++ )
          for ( j=0;j < mat->cols; j++ )
          {
            mat->array.N[i][j] = mat->kgv;
          }
      }
      else
      {
        for ( i=0; i < mat->rows; i++ )
          for ( j=0;j < mat->cols; j++ )
          {
            mat->array.N[i][j] = mat->kgv;
            Normal2( &mat->array.SZ[i][j], &mat->array.N[i][j]);
          }
      }
      result = 0;
    }
  }   
  return result;
}

/**************************************************************************\
@---------------------------------------------------------------------------
@ result = rat2kgv( matrix_TYP *mat );
@
| wandelt eine Matrix in Bruchdarstellung in eine 
| Matrix mit kgv-Darstellung um.
| "result" ist 0, falls alles glattging und nimmt den Wert
| -2, falls mat->prime != 0,
@ changes a matrix with mat->array.N != NULL to a matrix
@ with mat->kgv != 0 and mat->array.N = NULL
@ The result is 0, if it worked and -2, if mat->prime != 0
@---------------------------------------------------------------------------
@
\**************************************************************************/

int 
rat2kgv (matrix_TYP *mat)
{
int i, j, result;

  if ( mat->prime ) /* this would be quite idiotic */
  { 
    result = -2;
  } 
  else if ( mat->array.N == NULL )
  {
/*
 * really nothing to do
 */
    mat->flags.Integral = mat->kgv == 1;
    result = 0;
  } 
  else
  {     
/*
 *   kuerzen um Overflow Gefahr zu verringern
 */
    for ( i=0;i < mat->rows; i++)
      for ( j=0;j < mat->cols; j++)
      {                                
        Normal2( &mat->array.SZ[i][j], &mat->array.N[i][j] );
      }
/*
 *   kgv aller Nenner berechnen.
 */
    mat->kgv= 1;
    for ( i=0;i < mat->rows; i++)
      for ( j=0;j < mat->cols; j++)
      {                                
        mat->kgv = KGV( mat->kgv, mat->array.N[i][j] );
      }
/*
 *  alle Brueche erweitern
 */   
    if ( mat->kgv > 1 )
    {
      for ( i=0;i < mat->rows; i++)
        for ( j=0;j < mat->cols; j++)
        {
          mat->array.SZ[i][j] *= (mat->kgv/mat->array.N[i][j]);
        }                      
    }
    free2dim ( (char **)mat->array.N, mat->rows );
    mat->array.N = NULL;
    result = 0;
  }
  return result;
}

/*}}}  */
