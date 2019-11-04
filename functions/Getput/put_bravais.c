#include "typedef.h"
#include "getput.h"
#include "matrix.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: put_bravais.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ void fput_bravais(outfile, G, comment)
@ FILE *outfile;
@ bravais_TYP *G;
@ char *comment;
@  the same as put_bravais, but outfile must already have         |
@  been opened                                                    |
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
fput_bravais (FILE *outfile, bravais_TYP *G, const char *comment)
{
  int i;
  
/*---------------------------------------------------------------*\
| print header line, i.e. gen_no, form_no, zentr_no, normal_no    |
| and cen_no                                                      |
\*---------------------------------------------------------------*/

  if(G->form_no == 0 && G->zentr_no == 0 && G->normal_no == 0 && G->cen_no == 0)
    fprintf(outfile, "#%d  %% number of generators", G->gen_no);
  else
  {
     fprintf(outfile, "#g%d ", G->gen_no);
     if(G->form_no != 0)
       fprintf(outfile, "f%d ", G->form_no);
     if(G->zentr_no != 0)
       fprintf(outfile, "z%d ", G->zentr_no);
     if(G->normal_no != 0)
       fprintf(outfile, "n%d ", G->normal_no);
     if(G->cen_no != 0)
       fprintf(outfile, "c%d", G->cen_no);
  }
  if(comment != NULL)
  {
     fprintf(outfile, "%% %s", comment);
  }
  fprintf(outfile, "\n");

/*--------------------------------------------------------------------*\
|  print the matrices                                                  |
\*--------------------------------------------------------------------*/
  for(i=0; i<G->gen_no; i++)
  { 
#if 0
    fprintf( outfile, "putting generator %d\n",i);
#endif
    fput_mat(outfile, G->gen[i], "generator", 0);
#if 0
    fprintf( outfile, "done\n");
#endif
  }
  for(i=0; i<G->form_no; i++)
  {
    Check_mat(G->form[i]);
    fput_mat(outfile, G->form[i], "invariant form", 2);
  }
  for(i=0; i<G->zentr_no; i++)
    fput_mat(outfile, G->zentr[i], "zentr", 0);
  for(i=0; i<G->normal_no; i++)
  {
    Check_mat(G->normal[i]);
    fput_mat(outfile, G->normal[i], "generator of normalizer", 2);
  }
  for(i=0; i<G->cen_no; i++)
  {
    Check_mat(G->cen[i]);
    fput_mat(outfile, G->cen[i], "generator of centralizer", 2);
  }

/*--------------------------------------------------------------------*\
| print order of the bravais-group                                     |
\*--------------------------------------------------------------------*/
  fput_order(outfile, G->divisors, G->order);

}



/**************************************************************************\
@---------------------------------------------------------------------------
@ void put_bravais(G, filename, comment)
@ bravais_TYP *G;
@ char *filename;
@ char *comment;
@
@ prints the bravais_TYP G to the file with name 'filename'.
@ If 'filename' == NULL, G is printed to standard output.
@ comment is a string to write comments that will be ignored if
@ the output is used as input for another programm
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
put_bravais (bravais_TYP *G, const char *filename, const char *comment)
{
  int i;
  FILE *outfile;
  
/*------------------------------------------------------------*\
| Open output file												|
\*------------------------------------------------------------*/
if ( filename != NULL ) {
	if ( (outfile = fopen (filename, "w")) == NULL ) {
		fprintf (stderr, "put_bravais: Could not open %s\n", filename);
		exit (4);
		}
	}
else 
	outfile = stdout;


/*--------------------------------------------------------------------*\
| print header line                                                    |
\*--------------------------------------------------------------------*/

  /* changed by tilman on request to have an coherent bravais_TYP 
  if(G->form_no == 0 && G->zentr_no == 0 && G->normal_no == 0)
    fprintf(outfile, "#%d  %% number of generators", G->gen_no);
  else */
  { 
     fprintf(outfile, "#g%d ", G->gen_no);
     if(G->form_no != 0)
       fprintf(outfile, "f%d ", G->form_no);
     if(G->zentr_no != 0)
       fprintf(outfile, "z%d ", G->zentr_no);
     if(G->normal_no != 0)
       fprintf(outfile, "n%d ", G->normal_no);
     if(G->cen_no != 0)
       fprintf(outfile, "c%d ", G->cen_no);
  }
  if(comment != NULL)
  {
     fprintf(outfile, "%% %s", comment);
  }
  fprintf(outfile, "\n");
/*--------------------------------------------------------------------*\
|  print the matrices                                                  |
\*--------------------------------------------------------------------*/
  for(i=0; i<G->gen_no; i++)
    fput_mat(outfile, G->gen[i], "generator", 0);
  for(i=0; i<G->form_no; i++)
  {
    Check_mat(G->form[i]);
    fput_mat(outfile, G->form[i], "invariant form", 2);
  }
  for(i=0; i<G->zentr_no; i++)
    fput_mat(outfile, G->zentr[i], "zentr", 0);
  for(i=0; i<G->normal_no; i++)
  {
    Check_mat(G->normal[i]);
    fput_mat(outfile, G->normal[i], "generator of normalizer", 2);
  }
  for(i=0; i<G->cen_no; i++)
  {
    Check_mat(G->cen[i]);
    fput_mat(outfile, G->cen[i], "generator of centralizer", 2);
  }

/*--------------------------------------------------------------------*\
| print the order of the bravais_group                                 |
\*--------------------------------------------------------------------*/
  fput_order(outfile, G->divisors, G->order);

/*--------------------------------------------------------------------*\
|  close output-file                                                   |
\*--------------------------------------------------------------------*/
if ( outfile != stdout ) {
	fclose (outfile);
	}
}
/*{{{}}}*/
