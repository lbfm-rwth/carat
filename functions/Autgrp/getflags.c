#error THIS FILE IS NOT USED

#include"typedef.h"
/*****	This file contains some routines for input/output	*****/


/*****************************************************\
|	gets the options from the command line
\*****************************************************/
void 
getflags (flagstruct *fl, int *options)
{
	int	i;

/* depth for the scalar product combinations */
	fl->DEPTH = 0;
/* only the point stabilizer of the first STAB basis-vectors will be computed */
	fl->STAB = 0;
/* flag that Bacher-polynomials will be used */
	fl->BACH[0] = 0;
/* flag that the depth for the Bacher-polynomials is given as an argument, 
   default is 1 */
	fl->BACH[1] = 0;
/* depth for the Bacher-polynomials */
	fl->BACHDEP = 0;
/* flag that the scalar product for the Bacher-polynomials is given as an 
   argument, default is 1/2*norm of the vector */
	fl->BACH[2] = 0;
/* scalar product for the Bacher-polynomials */
	fl->BACHSCP = 0;
/* flag that generators will be read */
	fl->GEN = 0;
/* flag that every new generator is immediately written on the file 
   AUTO.tmp */
	fl->PRINT = 0;
/* flag that the vectors will be read instead of calculated by the program */
	fl->VEC = 0;
/* flag for the output-style: 0 means ASCII, 1 GAP-format, 2 MAGMA-format */
	fl->OUTPUT = 0;
/* scan through the arguments */
  for (i = 1; i < argc; ++i)
  {
    /* every option should start with a '-' */
    if ((str = strchr(argv[i], '-')) != NULL)
    {
      if (strlen(str) <= 1)
        fprintf(stderr, "unknown option %s: ignored\n", str);
      else if (str[1] == 'D')
      /* option -Dn where n is some non-negative integer for depth of
         scalar product  combinations */
      {
	if (strlen(str) <= 2  ||  strcspn(str+2, "0123456789") > 0)
	{
          fprintf(stderr, "Error: no non-negative integer specified with -D option\n");
	  exit (3);
	}
	else
		fl->DEPTH = atoi(str+2);
      }
      else if (str[1] == 'S')
      /* option -Sn where n is some non-negative integer for 
         n-point stabilizer */
      {
	if (strlen(str) <= 2  ||  strcspn(str+2, "0123456789") > 0)
	{
		fprintf(stderr, "Error: no non-negative integer specified with -S option\n");
		exit (3);
	}
	else
		fl->STAB = atoi(str+2);
      }
      else if (strlen(str) >= 3  &&  strncmp(str, "-BD", 3) == 0)
      /* option -BDn where n is some non-negative integer for depth of 
         Bacher-polynomials */
      {
	if (strlen(str) <= 3  ||  strcspn(str+3, "0123456789") > 0)
	{
         fprintf(stderr, "Error: no non-negative integer specified with -BD option\n");
         exit (3);
	}
	else
	{
		fl->BACHDEP = atoi(str+3);
		fl->BACH[0] = 1;
		fl->BACH[1] = 1;
	}
      }
      else if (strlen(str) >= 3  &&  strncmp(str, "-BS", 3) == 0)
      /* option -BSn where n is some integer for scalar product of 
         Bacher-polynomials */
      {
	if (strlen(str) <= 3  ||  strcspn(str+3, "-0123456789") > 0)
	{
          fprintf(stderr, "Error: no integer specified with -BS option\n");
          exit (3);
	}
	else
	{
		fl->BACHSCP = atoi(str+3);
		fl->BACH[0] = 1;
		fl->BACH[2] = 1;
	}
      }
      else if (strlen(str) == 2  &&  str[1] == 'B')
      /* option -B indicates that Bacher-polynomials will be used */
		fl->BACH[0] = 1;
      else if (strlen(str) == 2  &&  str[1] == 'G')
      /* option -G indicates that some generators can be read from the
         input stream */
	fl->GEN = 1;
      else if (strlen(str) == 3  &&  str[1] == 'V')
      /* option -V1 indicates that the short vectors for the first lattice
         are read  from the input stream, option -V2 that the short vectors
         for the second lattice are read and option -V3 that the short
         vectors for both lattices are read
         (-V2 only makes sense in the isometry-program) */
      {
	if (strpbrk(str+2, "123") == NULL)
	{
          fprintf(stderr, "Error: no integer between 1 and 3 specified with -V option\n");
	  exit (3);
	}
	else
		fl->VEC = atoi(str+2);
      }
      else if (strlen(str) == 2  &&  str[1] == 'V')
      /* option -V indicates that the short vectors are read from the
         input stream */
      {
	if (fl->VEC == 0)
		fl->VEC = 3;
      }
      else if (strlen(str) == 2  &&  str[1] == 'P')
      /* option -P indicates that new generators will be written
         to the file AUTO.tmp immediately */
	fl->PRINT = 1;
      else if (strlen(str) == 3  &&  str[1] == 'O')
      /* option -OG indicates that the output is converted to GAP-format,
         -OM that it is converted to MAGMA-format */
      {
	if (strpbrk(str+2, "GM") == NULL)
	{
          fprintf(stderr, "Error: can only convert to GAP (-OG) or MAGMA (-OM)\n");
          exit (3);
	}
	else if (str[2] == 'G')
		fl->OUTPUT = 1;
	else if (str[2] == 'M')
		fl->OUTPUT = 2;
      }
      else 
	fprintf(stderr, "unknown option %s: ignored\n", str);
     }
   }
/* if Bacher-polynomials are to be used and no depth is given it is set to
   the default-value 1 */
	if (fl->BACH[0] == 1  &&  fl->BACH[1] == 0)
		fl->BACHDEP = 1;
}
