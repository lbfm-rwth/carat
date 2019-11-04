#include "typedef.h"
#include "tools.h"
#include "matrix.h"
#include "getput.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: put_mat.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
/*{{{}}}*/
/*{{{  put_mat*/
/*
@
@ void put_mat( mat, file_name, comment, options );
@ void fput_mat (outfile, mat, comment, options)
@ 
@ result: void
@ arguments: 
@    matrix_TYP *mat  -- the matrix to be written
@    char * file_name -- the name of the file the matrix will be written to
@                        if file_name == NULL, we write to stdout
@    char * comment   -- optional comment written to the output-file
@  unsigned int options -- options to controll the output of the matrix
@
@     Currently, there are two options. The one is called PM_RATIONAL and
@     has a value of 0x00000001 (defined in ../../include/matrix.h ). It
@     causes the matrix to be written with numerator and denominator for
@     each entry if it isn't
@     integral. The flag is still ignored for matrices that have array.N
@     allocated. If PM_NORMALIZED is not set, even rational matrices are
@     written in lcd-representation]
@     The other flag is called PM_SHORTCUT and causes the flag-entries
@     of the matrix not to be ignored, i.e. the matrix is not printed as 
@     an NxM array, if it is symmetric, diagonal or even scalar.
@
@     This is the first time since I started to compet with this confused
@     source-code that I discover that someone used bit-flags, wich is a
@     nice thing, though the one forgot to #define and document them.
@
@     These flags should be used the following way:
@
@         put_mat( ... , PM_RATIONAL );
@ or      put_mat( ... , PM_SHORTCUT );             
@ or      put_mat( ... , PM_SHORTCUT | PM_RATIONAL );
@ or just put_mat( ... , 0UL );
@
@   Note: "0UL" means just "(unsigned long)0"
@-------------------------------------------------------------------------
@
 */
void 
put_mat (matrix_TYP *mat, const char file_name[], const char comment[], unsigned long options)
{
FILE *outfile;




  /*
   * Open output file
   */
  if ( file_name != NULL ) {
    if ( (outfile = fopen (file_name, "w")) == NULL ) {
      perror("put_mat: Error in fopen()");
      exit (4);
    }
  } else {
    outfile = stdout;
  }
  fput_mat (outfile, mat, comment, options);
  if ( outfile != stdout ) {
    fclose (outfile);
#if 0
    fprintf (stderr, "%s written to %s\n", comment, file_name);
#endif
  }
  fflush(stdout);
}

/*}}}  */
/*{{{  fput_mat*/


/*--------------------------------------------------------------------*\
|  the same as put_mat, but output_file must already have been         |
| opened                                                               |
\*--------------------------------------------------------------------*/

void 
fput_mat (FILE *outfile, matrix_TYP *mat, const char comment[], unsigned long options)
{  
int *max_lenZ;
flag_TYP flags;
rational d, *max_lenQ;
boolean Normalized;
int i, j, str_len;
char string[132];
char format[32];

  Normalized = !(options & PM_RATIONAL);
  /*
   * claus: I've got to correct that one day! Normalized means to
   * print out a matrix with kgv.
   */
  if ( mat->array.N != NULL ) Normalized = 0;
  
  flags = mat->flags;
  if ( !(options & PM_SHORTCUT ) )
  {
    flags.Symmetric =
    flags.Diagonal  =
    flags.Scalar    = FALSE;
  }
  
  /*
   * Print header line                        
   */
  fprintf (outfile, "%d", mat->rows);
  if ( mat->rows != mat->cols ) {
    fprintf (outfile, "x%d", mat->cols);
  } else {
    if ( flags.Symmetric ) {
      if ( flags.Diagonal ) {
        fprintf (outfile, "d");
        if ( flags.Scalar ) {
          fprintf (outfile, "0");
        } else { 
          fprintf (outfile, "1"); /* ! Scalar */
        }
      } else { 
        fprintf (outfile, "x0");  /* !Diagonal */
      }
    } /* Symmetric */
  }
  if ( !flags.Integral && mat->array.N == NULL ) {
    fprintf (outfile, "\t/");
    if ( Normalized ) {
      fprintf (outfile, "%d\t", mat->kgv);
    } else {
      fprintf (outfile, "0\t" );
    }      
  }
  fprintf (outfile, "\t%% %s\n", comment);
  
  /*
   * print matrix                         
   */
  if ( flags.Integral || Normalized ) {
    if ( flags.Diagonal ) {
      if ( flags.Scalar ) {
        fprintf (outfile, "%d\n", mat->array.SZ[0][0]);
      } else {
        for ( i  = 0 ; i < mat->rows; i++  ) {
          fprintf (outfile, "%d ", mat->array.SZ[i][i]);
        }
        fprintf (outfile, "\n");
      }      /* !Scalar   */
    } else { /* !Diagonal */
      max_lenZ = (int *)malloc((unsigned)mat->cols*sizeof(int));
      for ( i = 0 ; i < mat->cols; i++  ) {
        max_lenZ[i] = 0;
      }
      for ( i = 0; i < mat->rows; i++  ) {
        for(j=0; j<(flags.Symmetric ? i+1 : mat->cols); j++) {
          str_len=sprintf(string,"%d",mat->array.SZ[i][j]);
          if ( str_len >= max_lenZ[j] ) {
            max_lenZ[j] = str_len + 1;
          }
        }
      }
      for (i = 0; i < mat->rows; i++) {
        for (j=0; j < (flags.Symmetric ? i+1 : mat->cols); j++)
        fprintf(outfile,"%*d",max_lenZ[j],mat->array.SZ[i][j]);
        fprintf (outfile, "\n");
      }
      free (max_lenZ);
    }
  } else { /* !Integral && !Normalized */
    if ( flags.Diagonal ) {
      if ( flags.Scalar ) {
        d.z = mat->array.SZ[0][0];
        if ( mat->array.N != NULL ) {
          d.n = mat->array.N[0][0];
        } else {
          d.n = mat->kgv;
        }
        Normal (&d);
        fprintf (outfile, "%d/%d\n", d.z, d.n);
      } else {
        for ( i = 0 ; i < mat->cols; i++ ) {
          d.z = mat->array.SZ[i][i];
          if ( mat->array.N != NULL ) {
            d.n = mat->array.N[i][i];
          } else {
            d.n = mat->kgv;
          }
          Normal (&d);
          if ( d.n == 1 ) {
            fprintf (outfile, " %d", d.z);
          } else {
            fprintf (outfile, " %d/%d", d.z, d.n );
          }
        }
        fprintf (outfile, "\n");
      }
    } else { /* !Diagonal */
      max_lenQ=(rational *)malloc(mat->cols*sizeof(rational));
      for (i = 0; i < mat->cols; i++ ) {
        max_lenQ[i] = Zero;             
      }
      for (i = 0; i < mat->rows; i++) {
        for( j = 0; j < (flags.Symmetric ? i+1 : mat->cols); j++ ) {
          d.z = mat->array.SZ[i][j];
          if ( mat->array.N != NULL ) {
            d.n = mat->array.N[i][j];
          } else {
            d.n = mat->kgv;
          }
          Normal (&d);
          if((str_len=sprintf(string,"%d",d.z))>=max_lenQ[j].z) {
            max_lenQ[j].z = str_len + 1;                         
          }
          if((str_len=sprintf(string,"%d",d.n))>max_lenQ[j].n ) {
            max_lenQ[j].n = str_len;
          }
        }
      }
      for ( i = 0; i < mat->rows; i++) {
        for(j=0;j<(flags.Symmetric ? i+1 : mat->cols); j++ ) {
          d.z = mat->array.SZ[i][j];
          if ( mat->array.N != NULL ) {
            d.n = mat->array.N[i][j];
          } else {
            d.n = mat->kgv;
          }
          Normal (&d);
          if ( d.n == 1 ) {
            sprintf(format,"%%%dd",max_lenQ[j].z);
            fprintf(outfile,"%*d",max_lenQ[j].z+max_lenQ[j].n+1,d.z);
          } else {
            sprintf(format,"%%%dd/%%-%dd",max_lenQ[j].z,max_lenQ[j].n);
            fprintf(outfile,format,d.z,d.n);
          }                           
        }
        fprintf (outfile, "\n");
      }
      free (max_lenQ);
    } 
  }
}
/*}}}  */
