#include "typedef.h"
#include "autgrp.h"
#include "getput.h"
#include "symm.h"
#include "datei.h"
#include "matrix.h"

int SFLAG;
int INFO_LEVEL;

int main (int argc, char *argv[])
{

	matrix_TYP **F, **Erz, *SV;
        int i, Fanz, Erzanz, options[6];
        bravais_TYP *G;
        FILE *infile;
        int Fmax;

        read_header(argc, argv);

        if (is_option('h') && optionnumber('h')==12){
           SFLAG = 1;
        }

        if(FILEANZ != 1)
        {
           printf("Usage: %s 'file' [-d=n] [-s=n] [-b] [-b=n] [-B=n] [-v] [-g] [-p]\n",argv[0]);
           printf("\n");
           printf("file: matrix_TYP containing a set of integral NxN-matrices, the first\n");
           printf("      of which is symmetric and positive definite.\n");
           printf("\n");
           printf("Computes a generating set for the (finite) group of all g in GL_N(Z) with\n");
           printf("g^Tr * F * g = F for all F in file.\n");
           printf("NOTE: if all F are symmetric, this is a Bravais group, otherwise a\n");
           printf("      generalised Bravais group (relevant for Bravais groups in families\n");
           printf("      like 2-1',2-1' etc.).\n");
           printf("\n");
           printf("NOTE: CARAT was developed for crystallographic groups in dimensions\n");
           printf("      up to 6. Most algorithms also work in higher dimensions.\n");
           printf("      However, integer overflow is not trapped in general.\n");
           printf("\n");
           printf("Options:\n");
           printf("-d=n    : Depth up to which scalar products are calculated. The value\n");
           printf("          should be small. Usefull if the automorphism group is expected\n");
           printf("          to be small.\n");
           printf("-s=n    : The n-point stabilizer with respect to different basis will be\n");
           printf("          calculated.\n");
           printf("-b=n    : Use Bacher polynomials up to depth n.\n");
           printf("-B=n    : Use Bacher polynomials with vectors having scalar product n\n");
           printf("-v,-g   : Read additional data from 'file'. If -v is given the program\n");
           printf("          assumes that the short vectors of the first form in 'file'\n");
           printf("          are given below the forms.\n");
           printf("          If -g is given, the program assumes known generators for\n");
           printf("          the automorphism group to be given below any other information in\n");
           printf("          'file'.\n");
           printf("-p      : Write additional output to the file AUTO.tmp\n");
           printf("\n");
           printf("Cf. Short, Shortest\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }
       /*------------------------------------------------------------*\
       |  Open input file
       \*------------------------------------------------------------*/
       if ( (infile = fopen (FILENAMES[0], "r")) == NULL )
       {
          fprintf (stderr, "Could not open input-file %s\n", FILENAMES[0]);
          exit (4);
       }
       /*------------------------------------------------------------*\
       |  Read the input
       \*------------------------------------------------------------*/
        F = fmget_mat(infile, &Fanz);
        if(is_option('v') == TRUE)
          SV = fget_mat(infile);
        else
        {
          Fmax = F[0]->array.SZ[0][0];
          for(i=1;i<F[0]->cols;i++)
          {
             if(F[0]->array.SZ[i][i] > Fmax)
               Fmax = F[0]->array.SZ[i][i];
          }
          SV = short_vectors(F[0], Fmax, 0, 0, 0, &i);
        }
        if(is_option('g') == TRUE)
           Erz = fmget_mat(infile, &Erzanz);
        else
        {  Erz = NULL, Erzanz = 0;}
       /*------------------------------------------------------------*\
       | Close input file
       \*------------------------------------------------------------*/
       if ( infile != stdin )
         fclose (infile);

       /*------------------------------------------------------------*\
       | Read the options
       \*------------------------------------------------------------*/

/*-------------------------------------------------------------------*\
| options is a pointer to integer (of length 6)
| The possible options are encoded in the following way:
| options[0]:	The depth, up to wich scalar product combinations
|		shall be calculated. The value should be small.
|		options[0] > 0 should be used only, if the automorphismn
|		group is expected to be small (with respect to the number
|		of shortest vectors).
| options[1]:	The n-point stabiliser with respect to different basis
|		will be calculated.
| options[2]:	If options[2] = 1, additional output is written to the
|               file AUTO.tmp
| options[3]:   If options[3] = 1, Bacher polynomials are used. 
|		If options[3] = 2, Bacher polynomial are used up to a depth
|		                   specified in options[4].
|		If options[3] = 3, Bacher polynomials are used, using 
|		                   combinations of vectors having the scalar
|		                   product specified in options[5]
|		options[3] = 4 is the combination of options[3] = 2 and
|                              options[3] = 3.
| options[4]:	A natural number number  or zero (if options[3] = 2 or 4)
| options[5]:	An integral number (if options[3] = 3 or 4)
|
|	It is possible to use NULL for options,
|	in this case option is assumed to be [0,0,0,0,0,0]
\*************************************************************************/

   for(i=0;i<6;i++)
     options[i] = 0;
  if(is_option('d') == TRUE)
     options[0] = optionnumber('d');
  if(is_option('s') == TRUE)
     options[1] = optionnumber('s');
  if(is_option('p'))
     options[2] = 1;
  if(is_option('B')  || is_option('b'))
     options[3] = TRUE;
  if(optionnumber('b') > 0)
  {
     options[4] = optionnumber('b');
     options[3] = 2;
  }
  if(optionnumber('B') != 0)
  {
     options[5] = optionnumber('B');
     options[3] = 3;
  }
  if(optionnumber('B') != 0 && optionnumber('b') > 0)
     options[3] = 4;


 G = autgrp(F, Fanz, SV, Erz, Erzanz, options);

 put_bravais(G, NULL, "");


 free_bravais(G);
 free_mat(SV);
 for (i=0;i<Fanz;i++) free_mat(F[i]);
 free(F);
 for (i=0;i<Erzanz;i++) free_mat(Erz[i]);
 if (Erzanz != 0) free(Erz);

 if (is_option('h') && optionnumber('h')==12){
    pointer_statistics(0,0);
 } 
 return 0;
}
