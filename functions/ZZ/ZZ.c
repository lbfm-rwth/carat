#include <signal.h>
#include <stdlib.h>
#include "typedef.h"
#include "voronoi.h"
#include "matrix.h"
#include "tools.h"
#include "ZZ_P.h"
#include "ZZ_cen_fun_P.h"

boolean QUIET = FALSE, TEMPORAER = FALSE, SHORTLIST = FALSE, NURUMF = FALSE,
  U_option = TRUE, G_option = TRUE, LLLREDUCED = FALSE;

int ZCLASS;

int COUNTER = 0;
int NUMBER = 1000, SUBDIRECT = FALSE, LEVEL = 500;	/*Default-Werte der maximalen Anzahl der berechneten
							   Zentrierungen bzw. Anzahl der Level */
FILE *ZZ_temp, *ZZ_list;
int MAT_ALLOC = 0;
int constituents = 0;
int verbose = 0;

void ZZ_intern (Gram, data, tree)
     matrix_TYP *Gram;
     ZZ_data_t *data;
     ZZ_tree_t *tree;
{
    ZZ_node_t *act, *new;
    int i, j, k, n, d, di, end_num, act_anz;
    int ABBRUCH = FALSE;
    
    act = tree->root;
    do {
	for (i = 0; i < data->p_consts.k; i++) {
	    for (j = 0; j < data->p_consts.s[i]; j++) {
		d = ZZ_epimorphs (data, i, j);
		di = d;
		end_num = data->EnCo[i][j].hi + 1;
		if (d != 0) {
		    act_anz = 0;
		    for (k = 1; d > 0; d--, k *= end_num) {
			for (n = k; n < 2 * k; n++) {
			    act_anz++;
			    ZZ_pick_epi (data, n, i, j);
			    new = ZZ_center (data, act, i, j);
			    ABBRUCH = ZZ_ins_node (Gram, data, tree, 
						   act, new, i, j);
			}
		    }
		    act->anz_tg = act_anz;
		    
		    for (k = 0; k < di; k++) {
			free_mat (data->epi_base[k]);
		    }
		    free (data->epi_base);
		    data->epi_base = NULL;
		}
	    }
	}
    } while ((!ABBRUCH) && (ZZ_successor (data, &act)));
    ZZ_fput_data (data, tree, ABBRUCH);
}

void ZZ_transpose_array (int **array, int size)
{
	int i, j, ZZ_swap;

	for (i = 0; i < size; i++) {
		for (j = i + 1; j < size; j++) {
			ZZ_swap = array[i][j];
			array[i][j] = array[j][i];
			array[j][i] = ZZ_swap;
		}
	}
}

/******************************************************************************
@
@  (void)ZZ( group, gram, divisors, options );
@
@  or
@
@  (void)ZZ( group, gram, divisors, options, num_sublattices, dim1, dim2, ... )
@
@  bravais_TYP *group;
@  matrix_TYP  *gram;
@  int *divisors;
@  char *options;     
@  int num_sublattices, dim1, ...
@
@  group: bravais group, the algorithm applies the ZZ-algorithm with regard
@         to the prime-divisors of the order of the group.
@  gram:  a gram matrix of a quadratic form that is invariant under the group
@         This form must be totally anisotropic.
@  divisors: an array of integers that REALLY MUST HAVE 100 (ONEHUNDRET) 
@            entries. The index of the array specifies the prime number and the
@            contents of the array the multiplicity of the prime as divisor of
@            the group order. The ZZ function does not need the multiplicity
@            it just uses a subset of the prime divisors of the group order 
@            to compute irreducible constituents of the matrix representation 
@            of the group modulo those primes. Thus it is all right to set
@            the entries of the array to some non-zero number, e. g. if one 
@            wants to compute the irreducible constituents modulo 2 and 
@            modulo 5, then one has to write
@             
@              divisors[2] = 1; divisors[5] = 1;
@
@            All other entries must be set to zero. (at least those that belong 
@            to prime number indices). 
@
@            If one simply wants to compute the irreducible constituents for
@            all divisors of the group order, it suffices to specify the 
@            value NULL for the argument divisors.
@
@  options: a string of characters, (i.e. "n30l24rq"), supported are:
@  num_sublattices, dim1, ... : see option "p". Must only be specified if the
@                               option "p" was specified
@
@           "n#": terminate after computation of "#" centerings
@           "l#": terminate afer  "#" iterations of the algorithm
@           "r" : use the "lll"-algorithm on the Z-bases of the centerings
@           "q" : "quit-mode" - suppress messages on stdout/stderr
@           "v" : "verbose-mode" - opposite of "q"
@           "t" : create a file in the current directory that
@                 contains additional information. See also options "ugsb"
@                 This feature is on by default and the file name defaults
@                 to ZZ.tmp. If this option is given the argument
@                 "none", then no output file is created.
@           "u" : computes elementary divisors of the matrices that contain
@                 the bases of the invariant sublattices (i.e. of
@                 group->zentr[i]) and writes them to the
@                 outputfile.
@           "g" : computes the elementary divisors of the gram-matrices of
@                 the invariant sublattices and writes them to the output file.
@                 Implies option "t".
@           "s" : "shortlist" - prints in short form the messages that are
@                 suppressed by "q". Independent of "q" option
@           "b" : "print change of base only" -- suppresses the output of
@                 the gram-matrices and elementary devisors to the output file.
@           "p" : "projections" -- compute only thos sublattices, that have
@                 surjective projections on a given number of sublattices
@                 with a given dimension. The number of the sublattices is
@                 given by the first argument after "options", that is, by
@                 "num_sublattices". The dimension of the first sublattice
@                 is given by "dim1", the second by the argument "dim2" etc
@
@                 NOTE: the number of the "dim#" arguments MUST match the
@                       "num_sublattices"
@
@                 Given this information, the algorithm devides the original
@                 lattice in sublattices that are spanned by
@                 <e_1, ..., e_dim1>,
@                 <e_(dim1+1), ..., e_(dim1+dim2)>, etc
@   
@   The computed centerings are stored in "group->zentr" as Z-bases
@   (matrix_TYP **zentr). The number of centerings is stored in 
@   "group->zentr_no".
@
@   The functions exits by calling exit() in case of an error of if
@   the parameters are inconsistent.
@
@   NOTE: the bravais group is supposed to operate on a lattice of column
@         vectors. "SPALTENKONVENTION"
@
@  TODO:
@        write more documentation.
@
******************************************************************************/

static void scan_options (char *options, int *projections, FILE *outputfile)
{
	char *optp, *help;
	int num_proj;
	char *temp_name = "ZZ.tmp";

	SHORTLIST = FALSE;
	ZCLASS = 0;
	NURUMF = FALSE;	/* only write the matrices of change of base to the 
			 * file "ZZ.tmp", i.e. the Z-bases of the invariant 
			 * sublattices.
			 */
	/*
	 * Default-Werte der maximalen Anzahl der berechneten
	 * Zentrierungen bzw. Anzahl der Level
	 */
	NUMBER = 1000;
	LEVEL = 500;
	LLLREDUCED = FALSE; /* so what? */
	TEMPORAER = TRUE;  /* schreibt einige Information ueber die gefundenen
			     * Zentrierungen in die Datei "ZZ.tmp"
			     */
	ZZ_temp = NULL;
	QUIET = FALSE;      /* make less noise */
	SUBDIRECT = FALSE;  /* berechnet nur solche Untergitter deren 
			     * Projektionen auf bestimmte Teilgitter surjektiv
			     * sind
			     */
	U_option = TRUE;    /* berechnet Elementarteiler der Basen der 
			     * inv. TG 
			     */
	G_option = TRUE;   /* berechnet Elementarteiler der Grammatrizen */

	if ((optp = strchr (options, (int) 'l')) != NULL) {
		LEVEL = atoi (optp + 1);
	}
	if ((optp = strchr (options, (int) 'n')) != NULL) {
		NUMBER = atoi (optp + 1);
	}
	if ((optp = strchr (options, (int) 'r')) != NULL) {
		LLLREDUCED = TRUE;
	}
	if ((optp = strchr (options, (int) 'q')) != NULL) {
		QUIET = TRUE;
	}
	if ((optp = strchr (options, (int) 'v')) != NULL) {
		QUIET = FALSE;
	}
	if ((optp = strchr (options, (int) 'p')) != NULL) {
		SUBDIRECT = TRUE;
		projections[0] = atoi (optp + 1);
		if (projections[0] == 0) {
			fprintf (stderr, "\"p\" option requires number of sublattices and their dimensions to be given.\n");
			exit (3);
		} else if (projections[0] > 6) {
			fprintf (stderr, "Maximal dimension is 6\n");
			exit(3);
		} else {
#if DEBUG
			printf ("%d\n", projections[0]);
#endif
			help = optp + 1;
			for (num_proj = 1;
			     num_proj <= projections[0];
			     num_proj++){
				if ((help = strchr (help, '/')) != NULL) {
					help++;
					projections[num_proj] = atoi (help);
					if (projections[num_proj] == 0) {
						fprintf (stderr, "\"p\" option requires number of sublattices and their dimensions to be given.\n");
						exit (3);
					}
				} else {
					fprintf (stderr, "\"p\" option requires number of sublattices and their dimensions to be given.\n");
					exit (3);
				}
			}
 		}
	}
	if ((optp = strchr (options, (int) 't')) != NULL) {
		if (outputfile == NULL) {
			TEMPORAER = FALSE;
		} else {
			TEMPORAER = TRUE;
			ZZ_temp = outputfile;
		}
	}
	if ((optp = strchr (options, (int) 'u')) != NULL) {
		U_option = FALSE;
	}
	if ((optp = strchr (options, (int) 'z')) != NULL) {
		ZCLASS = 1;
	}
	if ((optp = strchr (options, (int) 'Z')) != NULL) {
		ZCLASS = 2;
	}
	if ((optp = strchr (options, (int) 'g')) != NULL) {
		G_option = FALSE;
	}
	if ((optp = strchr (options, (int) 's')) != NULL) {
		SHORTLIST = TRUE;
	}
	if ((optp = strchr (options, (int) 'b')) != NULL) {
		NURUMF = TRUE;
	}
	if (U_option || G_option || NURUMF) {
		TEMPORAER = TRUE;
	}
	/*
	 *  globale Abbruchvariable
	 */
	COUNTER = 0; /* Anzahl der Zentrierungen */
	if (TEMPORAER && ZZ_temp == NULL) {
		if ((ZZ_temp = fopen (temp_name, "w+")) == NULL) {
			fprintf (stderr, "ZZ: scan_arg: Error, could not open temporary file \n");
			TEMPORAER = FALSE;
		}
	}
}

void ZZ(group, gram, divisors, options, outputfile)
     bravais_TYP *group;
     matrix_TYP *gram;
     int *divisors;
     char *options;
     FILE *outputfile;
{
	ZZ_data_t data;
	ZZ_tree_t tree;
	int i, result = 0;
	matrix_TYP *sylv, *help;
	int projections[7];

	/* zuerst wird getestet, ob alle noetigen Daten vorhanden sind. */
	if (group == NULL) {
		printf("ZZ: Error: no bravais group specified.\n");
		exit(3);
	} 
	if (gram == NULL) {
		printf("ZZ: Error: no gram matrix specified.\n");
		exit(3);
	}
	sylv = dsylv(gram);
	for (i=0;i <sylv->rows;i++) {
		if (sylv->array.SZ[i][i] <= 0) {
			fput_mat( stderr, gram, "Gram", 0);
			free_mat(sylv);
			printf("ZZ: Error: gram matrix not positiv definite");
			exit(3);
		}
	}
	free_mat(sylv);
	if (divisors == NULL) {
		divisors = group->divisors;
	}
	result = 0;
	for (i = 1; i < 100 && result == 0; i++) {
		if (divisors[i] != 0) {
#if DEBUG
			printf("divisors[%d] = %d\n", i, divisors[i]);
#endif
			result = divisors[i];
		}
	}
	if (result == 0) {
		printf("ZZ: Error: prime divisors not specified.\n");
		exit(3);
	}
	if (group->gen_no == 0) {
		printf("ZZ: Error: bravais group with unknown number of generators?\n");
		exit(3);
	} 
	if (group->gen == NULL) {
		printf("ZZ: Error: bravais group without generators?\n");
		exit(3);
	}
	scan_options (options, projections, outputfile);
	for (i = 0; i < group->gen_no; i++) {
		ZZ_transpose_array(group->gen[i]->array.SZ, 
				   group->gen[i]->cols);
	}
	/* Now check whether the gram matrix is really invariant under the
	 * group
	 *
	 * NOTE: scal_pr() is a ROW scalar product.
	 
	for (i = 0; i < group->gen_no ; i++) {
		int check;
		
		help = scal_pr(group->gen[i], gram, TRUE);
		check = cmp_mat(help, gram);
		free_mat(help);
		if (check != 0) {
			fput_mat( stdout, gram, "form", 0);
			printf("i: %d\n", i);
			printf("The gram matrix isn't invariant under the group.\n");
                        exit(3);
		}
	} */

	ZZ_get_data (group, gram, divisors, &data, &tree, projections);
	ZZ_intern (gram, &data, &tree);
	ZZ_put_data (group, &data, &tree);
	for (i = 0; i < group->zentr_no; i++) {
		ZZ_transpose_array(group->zentr[i]->array.SZ,
				   group->zentr[i]->cols);
	}
	for (i = 0; i < group->gen_no; i++) {
		ZZ_transpose_array (group->gen[i]->array.SZ, 
				    group->gen[i]->cols);
	}
	ZZ_free_data (&data);
}
