#include "ZZ.h"

extern int INFO_LEVEL;
extern int SFLAG;
extern int IDEM_NO;

extern char *optarg;
extern int optind;
FILE *outputfile;

void ZZ_usage (progname)
     char *progname;
{

        IDEM_NO = 0;

	fprintf (stderr,
 "Usage: %s -b -g -h -l <#level> -m -n <#number> -p -q -r -s -t out -u bravais-file gram-file\n\n",
		 progname);
    fprintf (stderr,
"-b  : Print only the matrices of change of base and their inverse.\n");
    fprintf (stderr,
"-g  : Do not compute elementary divisors of the gram matrix.\n");
    fprintf (stderr,
"-h  : Print this help\n");
    fprintf (stderr,
"-l #: Stop after reaching level #level (default #=%d).\n", LEVEL);
    fprintf (stderr,
"-n #: Stop after computation of #number \"Zentrierungen\" (default #=%d).\n", NUMBER);
    fprintf (stderr,
"-p<d0>/<d1>/<d2> ... : treat the lattice as a direct sum of <d0> sublattices\n");
    fprintf (stderr,
"      of dimensions <d1>, <d2> etc. (0 <= d0 <= 6) and compute only those\n");
    fprintf (stderr,
"      centerings that have surjective projections on them.\n");
    fprintf (stderr,
"-q  : Quiet mode. Suppress any messages to stdin/stdout.\n");
    fprintf (stderr,
"-r  : With ZZ_lll-reduction.\n");
    fprintf (stderr,
"-s  : Print less information.\n");
    fprintf (stderr, 
"-t  : Create an output file with additional information. The name of the\n");
    fprintf (stderr, 
"      output file defaults to stdout. Specifying \"none\" disables\n");
    fprintf (stderr,
"      writing to the output file\n");
    fprintf (stderr,
"-m  : Used for debugging, do not use!\n");
    fprintf (stderr,
"-u  : Do not compute elementary divisors of the change of base\n");
    fprintf (stderr, 
"\"bravais-file\" is a file containing a bravais group\n");
    fprintf (stderr,
"\"gram-file\" contains a positive definite invariant gram matrix.\n");
    fprintf (stderr, "\n");
}

int parse_options(argc, argv, options)
     int argc;
     char *argv[];
     char *options;
{
	int c;
	int errflag = 0;
	int i;
	char *p;

	*options++ = 't';
	outputfile = stdout;
	/*  Ok, it seems to be carat convention to allow a '=' sign
	 *  between the option and the argument.
	 *  We cope with this by simply replacing any occurences of a '=' sign 
	 *  with spaces.
	 */
	for (i=1; i < argc; i++) {
		for (p=argv[i]; *p != '\0'; p++) {
			if (*p == '=') {
				strcpy(p, p+1);
				break;
			}
		}
	}
	while ((c = getopt(argc, argv, "bghl:n:p:zZD:qrst:mu")) != -1) {
		switch(c) {
		case 'b':
		case 'g':
		case 'q':
		case 'r':
		case 's':
		case 'm':
		case 'u':
		case 'z':
		case 'Z':
			*options++ = (char)c;
			break;
		case 'l':
		case 'D':
		case 'n':
		case 'p':
			*options++ = (char)c;
			strcpy(options, optarg);
			options += strlen(optarg);
			break;			
		case 't':
			*options++ = (char)c;
			if (strcmp("none", optarg)) {
				if (strcmp(optarg, "-")) {
					outputfile = fopen(optarg, "w+");
				}
				if (outputfile == NULL) {
					perror("ZZprog: Error, could not open temporary file\n");
					errflag = 1;
				}
			} else {
				outputfile = NULL;
			}
			break;			
		case 'h':
		case '?':
		default:
			fprintf(stderr, "Unrecognized option\n");
			errflag++;
			break;
		}
	}
	return errflag == 0;
}

int main (argc, argv)
     int argc;
     char *argv[];
{
    bravais_TYP *group;
    matrix_TYP *gram = NULL;
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

	    for ( i = optind; i < argc; i++) {
		    printf("argv[%d] = %s\n", i, argv[i]);
	    }
    }
#endif
    switch (argc - optind) {
    case 1:
	    group = get_bravais (argv[optind]);
	    gram = group->form[0];
	    break;
    case 2:
	    group = get_bravais (argv[optind]);
	    gram = get_mat(argv[optind+1]);
	    break;
    default:
	    ZZ_usage (prog_name);
	    exit(0);
    }


#if DEBUG
    fput_mat( stderr, gram, "Form", 0);
#endif

    ZZ (group, gram, group->divisors, options, outputfile);

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
