#include "typedef.h"
#include "tools.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE:  malloc2dim.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/

/*
 |
 | tools/malloc2dim.c -- malloc & co fuer 2-dim. Arrays
 | exportiert die Funktionen
 |
 | calloc2dim
 | malloc2dim
 | memcpy2dim
 | memset2dim
 | free2dim
 |
*/

/*{{{}}}*/
/*{{{  calloc2dim*/
/*
@
@-------------------------------------------------------------------------
@ char **calloc2dim ( r, c, size )
@ int r,c,size;
|-- allokiert ein 2-dimesionales array
@
@ allocates a 2-dimensional array with 'r' rows and 'c' columns
@ the size of the entries in bytes is given by the argument 'size'.
| int r, int c: ZeilenxSpalten des Arrays
| int size: Groesse eines Eintrags in bytes
@
*/
char **calloc2dim(r,c,size)
int r,c,size;
{
char **new; 
int i, j;

  if( r == 0 || c == 0 || size == 0 )
    new= NULL;
  else 
  {
    new= (char **)malloc( r*sizeof(char *));
    if ( new != NULL )
    {
      i= 0;
      do
      {
        new[i]= (char *)calloc( c, size );
      } while ( new[i++] && i < r );
      if ( i != r )
      { 
        for ( j=0; j < i;j++) free ( (int *)new[i] );
        free ( (int *)new );
        new= NULL;
      }                         
    }
  }
  return new;
}

/*}}}  */
/*{{{  malloc2dim*/
/* 
@-------------------------------------------------------------------------
@ char **malloc2dim ( r, c, size )
@ int r,c,size;
| -- allokiert ein 2-dimesionales array
@
| int r, int c: ZeilenxSpalten des Arrays
| int size: Groesse eines Eintrags in bytes
@ allocates a 2-dimensional array with 'r' rows and 'c' columns
@ the size of the entries in bytes is given by the argument 'size'.
@
@-------------------------------------------------------------------------
*/
char **malloc2dim(r,c,size)
int r,c,size;
{ 
char **new; 
int i, j;

  if( r == 0 || c == 0 || size == 0 )
    new= NULL;
  else 
  {
    new= (char **)malloc( r*sizeof(char *));
    if ( new != NULL )
    {
      i= 0;
      do
      {
        new[i]= (char *)malloc( c*size);
      } while ( new[i++] && i < r );
      if ( i != r )
      { 
        for ( j=0; j < i;j++) free ( (int *)new[i] );
        free ( (int *)new );
        new= NULL;
      }                         
    }
  }
  return new;
}

/*}}}  */
/*{{{  memcpy2dim*/
/*
@-------------------------------------------------------------------------
@ void memcpy2dim ( dest, src, r, c, size )
@        copies a 2-dimesional array
@
@ char **dest: destination
@ char **src:  source
@ int r, int c: rows x columns of the arrays
@ int size: size of an entry in bytes
@
@-------------------------------------------------------------------------
*/
void memcpy2dim(dest, src, r,c,size)
char **dest, **src;
int r,c,size;
{ 
int i, j;

  if( !( r == 0 || c == 0 || size == 0 || dest == NULL || src == NULL ) )
  { 
    for ( i=0; i < r; i ++ ) memcpy(dest[i],src[i],c*size);
  }
}

/*}}}  */
/*{{{  memset2dim*/
/*
@-------------------------------------------------------------------------
@ void memset2dim ( dest, r, c, size, value )
@  initializes a 2-dimensional array
@
@ char **dest: destination
@ int r, int c: rows x columns of the arrays
@ int size: size of an entry in bytes
| char *value: Pointer auf den Eintrag, mit dem das Feld initialisiert 
|              werden soll
@ char *value: pointer to the entry, the array shall be initialized with
@
@-------------------------------------------------------------------------
*/
void memset2dim(dest, r, c, size, value)
char **dest;
int r,c,size;            
char *value;
{ 
int i, j;

  if( !( r == 0 || c == 0 || size == 0 || dest == NULL || value == NULL) )
  { 
    for ( j=0; j < c; j++ ) {
      memcpy( &dest[0][j*size], value, size );
    }
    for ( i=1; i < r; i ++ ) {
      memcpy( dest[i], dest[0], size*c );
    }
  }
}

/*}}}  */
/*{{{  free2dim*/
/*
@-------------------------------------------------------------------------
@ void free2dim ( old, rows )
| -- gibt ein 2-dimesionales array frei
@ frees a 2-dimensional array
@
@ char **old: pointer to the array
| int rows: Anzahl der Zeilen (obere Grenze des ersten Arrayindexes)
@ int rows: number of rows (upper bound for the first index of the array)
@
@-------------------------------------------------------------------------
*/
void free2dim(old, rows)
char **old;        
int rows;
{ 
int i;

  for (i=0; i < rows;i++)
  { 
    free( (int *)old[i] );
  }               
  free( (int *)old );
}

/*}}}  */

