#include "typedef.h"
#include "matrix.h"
#include "getput.h"


int main (int argc, char *argv[])
{

	matrix_TYP **F1, **F2, **Erz, *SV1, *SV2, *Iso, *tmp;
        int i, Fanz, F2anz, Erzanz, options[6];
        FILE *infile;
        int Fmax;

        extern char **FILENAMES;
        extern int FILEANZ;

	extern matrix_TYP *fget_mat ();
	extern matrix_TYP **fmget_mat ();
	extern matrix_TYP *short_vectors();
	extern matrix_TYP *isometry();

        read_header(argc, argv);
        if(FILEANZ != 2)
        {
           printf("Usage:  %s 'file1' 'file2 [-d=n] [-s=n] [-b] [-b=n] [-B=n] [-g] [-v] [-p]\n",argv[0]);
           printf("\n");
           printf("file1: matrix_TYP containing a set of m NxN-matrices, the first of which,\n");
           printf("       must be symmetric and positive definite.\n");
           printf("file2: matrix_TYP containing a set of m NxN-matrices, the first of which,\n");
           printf("       must be symmetric and positive definite.\n");
           printf("\n");
           printf("Computes a single g in GL_n(Z) with the property g^Tr * F_i_1 * g = F_i_2 for all\n");
           printf("F_i_1 in 'file1' and F_i_2 in 'file2'.\n");
           printf("\n");
           printf("Options:\n");
           printf("-d=n    : Depth up to which scalar products are calculated. The value\n");
           printf("          should be small.\n");
           printf("-s=n    : The n-point stabilizer with respect to different basis will be\n");
           printf("          calculated.\n");
           printf("-b=n    : Use Bacher polynomials up to deepth n.\n");
           printf("-B=n    : Use Bacher polynomials with vectors having scalar product n\n");
           printf("-v,-g   : Read additional data from 'file'. If -v is given the program\n");
           printf("          assumes that the short vectors of the first form in 'file'\n");
           printf("          are given below the forms.\n");
           printf("          If -g is given, the program assumes known generators for\n");
           printf("          the automorphism group to be given below any other information in\n");
           printf("          'file'.\n");
           printf("-p      : Write additional output to the file ISOM.tmp\n");
           printf("\n");
           printf("Cf. Short, Shortest, Aut_grp\n");
           if (is_option('h')){
              exit(0);
           }
           else{
              exit(31);
           }
        }
       /*------------------------------------------------------------*\
       |  Open input file1
       \*------------------------------------------------------------*/
       if ( (infile = fopen (FILENAMES[0], "r")) == NULL )
       {
          fprintf (stderr, "Could not open input-file %s\n", FILENAMES[0]);
          exit (4);
       }
       /*------------------------------------------------------------*\
       |  Read the input
       \*------------------------------------------------------------*/
        F1 = fmget_mat(infile, &Fanz);
        if(is_option('v') == TRUE)
          SV1 = fget_mat(infile);
        else
        {
          Fmax = F1[0]->array.SZ[0][0];
          for(i=1;i<F1[0]->cols;i++)
          {
             if(F1[0]->array.SZ[i][i] > Fmax)
               Fmax = F1[0]->array.SZ[i][i];
          }
          SV1 = short_vectors(F1[0], Fmax, 0, 0, 0, &i);

        }
       /*------------------------------------------------------------*\
       | Close input file
       \*------------------------------------------------------------*/
       if ( infile != stdin )
         fclose (infile);

       /*------------------------------------------------------------*\
       |  Open input file2
       \*------------------------------------------------------------*/
       if ( (infile = fopen (FILENAMES[1], "r")) == NULL )
       {
          fprintf (stderr, "Could not open input-file %s\n", FILENAMES[1]);
          exit (4);
       }
       /*------------------------------------------------------------*\
       |  Read the input
       \*------------------------------------------------------------*/
        F2 = fmget_mat(infile, &F2anz);
        if(F2anz != Fanz)
        {
          printf("Error: different number of matrices in 'isometry'\n");
          exit(3);
        }
        if(is_option('v') == TRUE)
          SV2 = fget_mat(infile);
        else
          SV2 = short_vectors(F2[0], Fmax, 0, 0, 0, &i);
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
|		If options[3] = 2, Bacher polynomial are used up to a deepth
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
  if(optionnumber('B') > 0)
  {
     options[5] = optionnumber('4');
     options[3] = 3;
  }
  if(optionnumber('B') > 0 && optionnumber('b') > 0)
     options[3] = 4;


 Iso = isometry(F1, F2, Fanz, SV1, SV2, Erz, Erzanz, options);
 if(Iso != NULL) {
   tmp = tr_pose(Iso);
   put_mat(tmp, NULL, "isometry", 0);
   free_mat(tmp);
 }
 else
   printf("The forms are not isometric\n");


 exit(0);
}
