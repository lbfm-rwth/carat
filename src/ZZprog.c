#include "ZZ.h"
#include "bravais.h"


static FILE *outputfile;

void 
ZZ_usage (char *progname)
{
     printf("Usage: %s 'file1'  ['file2'] [-b] [-g] [-h] [-l=<#level>] \n"
                                  ,progname);
     printf("       [-m] [-n=<#number>] [-p] [-q] [-r] [-s] [-t=out] [-u]\n");
     printf("\n");
     printf("file1: bravais_TYP containing a finite unimodular group G,\n");
     printf("       ending with the order of G.\n");
     printf("file2: (OPTIONAL) matrix_TYP containing the Gram matrix of a symmetric \n");
     printf("       positive definite G-invariant bilinear form. This form is\n");
     printf("       used for reduction purposes only. If this file is not\n");
     printf("       given, the program computes such a form. In particular\n");
     printf("       the forms possibly given in 'file1' are ignored.\n");
     printf("\n");
     printf("Calculates the the G-sublattices of the natural lattice Z^n of\n");
     printf("finite index, the prime divisors of which divide the order of G\n");
     printf("as given in 'file1'. Sublattices of proper multiples of Z^n are\n");
     printf("ignored.\n");
     printf("\n");
     printf("Options:\n");
     printf("-b      : Print only the matrices of change of base and their inverse.\n");
     printf("-g      : Do not compute elementary divisors of the gram matrix.\n");
     printf("-l=#    : Stop after reaching level #level (default #=500).\n");
     printf("-n=#    : Stop after computation of #number \"sublattices\" (default #=1000).\n");
     printf("-q      : Quiet mode. Suppress any messages to stderr.\n");
     printf("-r      : With LLL-reduction for the bases, cf. 'file2'.\n");
     printf("-s      : Print less information.\n");
     printf("-t='out': Create an output file with additional information. The name \n");
     printf("          of the output file defaults to stdout. Specifying \"none\" \n");
     printf("          disables writing to the output file\n");
     printf("-u      : Do not compute elementary divisors of the basis\n");
     printf("          transformations.\n");
     printf("\n");
     printf("Options for experts:\n");
     printf("-p<N>/<d1>/<d2>...<dN> : treat the lattice as a direct sum of <N>  \n");
     printf("                         sublattices of dimensions <d1>, <d2> etc.  \n");
     printf("                         (1 <= N, di <= 6) and compute only those\n");
     printf("                         sublattices that have surjective projections  \n");
     printf("                         onto each of the N component lattices.\n");
     printf("\n");
     printf("Cf. Order, QtoZ, Z_equiv, Q_equiv\n");

}

int parse_options(int argc,
                  char *argv[],
                  char *option)
{

   int i,
       errflag = 0;

   *option++ = 't';
   outputfile = stdout;

   if (is_option('b')){
      *option='b';
      option++;
   }
   if (is_option('g')){
      *option='g';
      option++;
   }
   if (is_option('q')){
      *option='q';
      option++;
   }
   if (is_option('r')){
      *option='r';
      option++;
   }
   if (is_option('s')){
      *option='s';
      option++;
   }
   if (is_option('m')){
      *option='m';
      option++;
   }
   if (is_option('u')){
      *option='u';
      option++;
   }
   if (is_option('z')){
      *option='z';
      option++;
   }
   if (is_option('Z')){
      *option='Z';
      option++;
   }

   if (is_option('l')){
      *option='l';
      option++;
      sprintf(option,"%d ",optionnumber('l'));
      option=strchr(option,' ');
   }
   if (is_option('n')){
      *option='n';
      option++;
      sprintf(option,"%d ",optionnumber('n'));
      option=strchr(option,' ');
   }

   if (is_option('p')){
      *option='p';
      option++;
      for (i=1;i<argc && strstr(argv[i],"-p") != argv[i];i++);
      sprintf(option,"%s ",argv[i]+2);
      option=strchr(option,' ');
   }
   if (is_option('t')){
      *option='t';
      option++;
      for (i=1;i<argc && strstr(argv[i],"-t") != argv[i];i++);
      outputfile = fopen(argv[i]+3, "w+");
      if (outputfile == NULL) {
         fprintf(stderr,"ZZprog: Error, could not open temporary file %s\n",
                       argv[i]+3);
         errflag = 1;
      }
   }

   return errflag == 0;
}

int 
main (int argc, char *argv[])
{
    bravais_TYP *group;
    matrix_TYP *gram = NULL,
               *ID;
    char *prog_name = argv[0];
    char options[256];

    /*  for somewhat obscure carat option handling. Why not use the well
     *  ANSI standard, i.e. getopt() :-(
     */
    read_header(argc, argv);

    if (is_option('D') && optionnumber('D') == 8){
       SFLAG = 1;
       INFO_LEVEL = optionnumber('D');
    }

    if (!parse_options(argc, argv, options)) {
	    ZZ_usage (prog_name);
	    exit(31);
    }

#if DEBUG
    printf("Options: %s\n", options);
    {
	    int i;

	    for ( i = 0; i < argc; i++) {
		    printf("argv[%d] = %s\n", i, argv[i]);
	    }
    }
#endif
    switch (FILEANZ) {
    case 1:
	    group = get_bravais (FILENAMES[0]);
            ID = init_mat(group->dim,group->dim,"i1");
	    gram = rform(group->gen,group->gen_no,ID,101);
            free_mat(ID);
	    break;
    case 2:
	    group = get_bravais (FILENAMES[0]);
	    gram = get_mat(FILENAMES[1]);
	    break;
    default:
	    ZZ_usage (prog_name);
	    exit(0);
    }


#if DEBUG
    fput_mat( stderr, gram, "Form", 0);
#endif

    ZZ (group, gram, group->divisors, NULL, options, outputfile, 0, 0);

#if DEBUG
    fprintf (stderr, "num_zentr: %d\n", group->zentr_no);
#endif
    /*  the following is not exactly needed as we do an exit() anyways.
     *  It's useful, however, for debugging.
     */
    cleanup_prime();
    free_bravais(group);
    if (gram != NULL) free_mat(gram);

    if (INFO_LEVEL == 8){
       pointer_statistics(0,0);
    }

    return 0;
}
