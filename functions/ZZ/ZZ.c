#include <signal.h>
#include <stdlib.h>
#include "typedef.h"
#include "voronoi.h"
#include "matrix.h"
#include "tools.h"
#include "ZZ_P.h"
#include "ZZ_cen_fun_P.h"
#include "ZZ_zclass_P.h"
#include "ZZ.h"
#include "longtools.h"

boolean QUIET = FALSE, TEMPORAER = FALSE, SHORTLIST = FALSE, NURUMF = FALSE,
  U_option = TRUE, G_option = TRUE, LLLREDUCED = FALSE;

int ZCLASS;
extern int *SUB_VEC;
extern matrix_TYP **PrI;

int COUNTER = 0;
int NUMBER = 1000, SUBDIRECT = FALSE, LEVEL = 500;	/*Default-Werte der maximalen Anzahl der berechneten
							   Zentrierungen bzw. Anzahl der Level */
FILE *ZZ_temp;
int MAT_ALLOC = 0;
int constituents = 0;
int verbose = 0;
ZZ_super_TYP **SUPER_info, *SUPER_INFO;

int IDEM_NO;


/*------------------------------------------------------------------------------- */
static int suche_mat(matrix_TYP *mat,
                     matrix_TYP **liste,
                     int anz)
{
   int i;

   for (i = 0; i < anz; i++){
      if (cmp_mat(mat, liste[i]) == 0)
         return(i);
   }
   return(-1);
}


/*------------------------------------------------------------------------------- */
void 
ZZ_intern (matrix_TYP *Gram, ZZ_data_t *data, ZZ_tree_t *tree, QtoZ_TYP *inzidenz)
{
    ZZ_node_t *act, *newnode, *nnn;
    int g, i, j, k, l, m, n, d, di, end_num, act_anz, flag = 0, nr, NEU, zahl, flagge;
    int ABBRUCH = FALSE;
    matrix_TYP *gitter, *X, *Li;


    act = tree->root;
    do {
       if (ZCLASS == 1 && GRAPH){
          flag = suche_mat(act->U, inzidenz->gitter, inzidenz->anz);
          if (flag == -1 && inzidenz->anz != 0){
             /* Das aktuelle Gitter kommt nicht vor in der Liste */
             fprintf(stderr,"ERROR 1 in ZZ_intern\n");
             exit(2);
          }
          if (flag == -1 && inzidenz->anz == 0){
             /* insert start-lattice */
             inzidenz->anz++;
             inzidenz->gitter[0] = copy_mat(act->U);
             inzidenz->tr_gitter[0] = tr_pose(inzidenz->gitter[0]);
             inzidenz->inv_tr_gitter[0] = mat_inv(inzidenz->tr_gitter[0]);
             inzidenz->entry[0] = (QtoZ_entry_TYP *)calloc(1024, sizeof(QtoZ_entry_TYP));
             flag = 0;
          }
       }
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
			    newnode = ZZ_center (data, act, i, j);
			    nr = 0;
			    NEU = 0;
			
			    if (ZCLASS == 1 && GRAPH){
 			       gitter = tr_pose(newnode->U);
 			    } else
 			       gitter = NULL;
 			    flagge = 0;
			    ABBRUCH = ZZ_ins_node (Gram, data, tree,
						   act, newnode, i, j,
						   inzidenz, &nr, &NEU, &flagge, &g, &nnn);
		            if (ZCLASS == 1 && GRAPH){
			       if (NEU != 1){
			          /* not a new lattice for the graph or lattice in another zoo */
			          if (nr == -1){
			             /* lattice possibly in another zoo */
			             zahl = inzidenz->entry[flag][0].anz;
			             if (zahl == 0){
                                        inzidenz->entry[flag][0].I = (int *)calloc(1024,
                                                                           sizeof(int));
                                        inzidenz->entry[flag][0].J = (int *)calloc(1024,
                                                                           sizeof(int));
                                        inzidenz->entry[flag][0].flag = (int *)calloc(1024,
                                                                           sizeof(int));
                                        inzidenz->entry[flag][0].lattice = (matrix_TYP
                                                       **)calloc(1024, sizeof(matrix_TYP *));
                                        inzidenz->zoogitter[flag] = (matrix_TYP **)calloc(1024,
                                                                     sizeof(matrix_TYP *));			
			             }
			             inzidenz->entry[flag][0].I[zahl] = i;
                                     inzidenz->entry[flag][0].J[zahl] = j;
                                     inzidenz->entry[flag][0].flag[zahl] = -10;
                                     inzidenz->zoogitter[flag][zahl] = gitter;

                                     /* This isn't the correct conjugating matrix yet! */
                                     inzidenz->entry[flag][0].lattice[zahl] =
                                        copy_mat(inzidenz->inv_tr_gitter[flag]);
                                     gitter = NULL;
                                     inzidenz->entry[flag][0].anz++;
			          }
  			          else{
			             /* not a new lattice */
			             nr++;
			             zahl = inzidenz->entry[flag][nr].anz;
			             if (zahl == 0){
                                        inzidenz->entry[flag][nr].I = (int *)calloc(1024,
                                                                           sizeof(int));
                                        inzidenz->entry[flag][nr].J = (int *)calloc(1024,
                                                                           sizeof(int));			
                                        inzidenz->entry[flag][nr].flag = (int *)calloc(1024,
                                                                           sizeof(int));			
                                        inzidenz->entry[flag][nr].lattice = (matrix_TYP
                                                       **)calloc(1024, sizeof(matrix_TYP *));
                                     }
 			             inzidenz->entry[flag][nr].I[zahl] = i;
                                     inzidenz->entry[flag][nr].J[zahl] = j;
                                     inzidenz->entry[flag][nr].flag[zahl] = flagge;
                                     inzidenz->entry[flag][nr].anz++;

                                     switch (flagge){
                                        case 1:
/*
                                           inzidenz->entry[flag][nr].lattice[zahl] =
                                                mat_mul(inzidenz->inv_tr_gitter[flag],
                                                        inzidenz->tr_gitter[nr - 1]);
fprintf(stderr, "F A L L   1!!!!!!!!!!!!!\n");
                                           break;
*/
                                        case 100:
                                           inzidenz->entry[flag][nr].lattice[zahl] =
                                                mat_mul(inzidenz->inv_tr_gitter[flag],
                                                        inzidenz->tr_gitter[nr - 1]);
                                           for (l = 0; l < gitter->rows; l++){
                                              for (m = 0; m < gitter->cols; m++){
                                                 inzidenz->entry[flag][nr].lattice[zahl]->array.SZ[l][m]
                                                    *= g;
                                              }
                                           }
                                           break;

                                        case 10:
                                           Li = mat_inv(gitter);
                                           X = konjugierende(Li, tree->root->col_group, nnn);
                                           free_mat(Li);
                                           if (X == NULL){
                                              fprintf(stderr, "ERROR 2 in ZZ_intern!\n");
                                              exit(9);
                                           }
                                           inzidenz->entry[flag][nr].lattice[zahl] =
                                                mat_mul(inzidenz->inv_tr_gitter[flag], gitter);
                                           mat_muleq(inzidenz->entry[flag][nr].lattice[zahl], X);
                                           free_mat(X);
                                           break;
                                     }
			          }
			       }
			       else{
			          /* new lattice in the graph */
                                  inzidenz->gitter[inzidenz->anz] = copy_mat(newnode->U);
                                  inzidenz->tr_gitter[inzidenz->anz] = tr_pose(newnode->U);
                                  inzidenz->inv_tr_gitter[inzidenz->anz] =
                                       mat_inv(inzidenz->tr_gitter[inzidenz->anz]);
                                  inzidenz->entry[inzidenz->anz] = (QtoZ_entry_TYP *)calloc(1024,
                                                                   sizeof(QtoZ_entry_TYP));
                                  inzidenz->anz++;
                                  inzidenz->entry[flag][inzidenz->anz].anz = 1;
                                  inzidenz->entry[flag][inzidenz->anz].I = (int *)calloc(1024,
                                                                           sizeof(int));
                                  inzidenz->entry[flag][inzidenz->anz].J = (int *)calloc(1024,
                                                                           sizeof(int));
                                  inzidenz->entry[flag][inzidenz->anz].flag = (int *)calloc(1024,
                                                                           sizeof(int));
                                  inzidenz->entry[flag][inzidenz->anz].lattice = (matrix_TYP
                                                       **)calloc(1024, sizeof(matrix_TYP *));
                                  inzidenz->entry[flag][inzidenz->anz].I[0] = i;
                                  inzidenz->entry[flag][inzidenz->anz].J[0] = j;
                                  inzidenz->entry[flag][inzidenz->anz].flag[0] = 0;
                                  inzidenz->entry[flag][inzidenz->anz].lattice[0] =
                                          mat_mul(inzidenz->inv_tr_gitter[flag],
                                                  inzidenz->tr_gitter[inzidenz->anz - 1]);
                               }
                               if (gitter != NULL)
                                  free_mat(gitter);
                            }
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
    nnn = NULL;
}




/******************************************************************************/
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

static void scan_options (const char *options, int *projections, FILE *outputfile)
{
	const char *optp, *help;
	int num_proj;
	const char *temp_name = "ZZ.tmp";

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
		} else if (projections[0] > 8) {
			fprintf (stderr, "Maximal dimension is 8\n");
			exit(3);
		} else {
#ifdef DEBUG
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

void
ZZ (bravais_TYP *group, matrix_TYP *gram, int *divisors, QtoZ_TYP *inzidenz, const char *options, FILE *outputfile, int super_nr, int konst_flag)
{
	ZZ_data_t *data;
	ZZ_tree_t *tree;
	int i, result = 0;
	matrix_TYP *sylv;
	int projections[9];
	ZZ_node_t *n, *t;
	ZZ_super_TYP *DATEN = 0;

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
#ifdef DEBUG
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
		
		matrix_TYP * help = scal_pr(group->gen[i], gram, TRUE);
		check = cmp_mat(help, gram);
		free_mat(help);
		if (check != 0) {
			fput_mat( stdout, gram, "form", 0);
			printf("i: %d\n", i);
			printf("The gram matrix isn't invariant under the group.\n");
                        exit(3);
		}
	} */
	
	
        data = (ZZ_data_t *)calloc(1, sizeof(ZZ_data_t));
	tree = (ZZ_tree_t *)calloc(1, sizeof(ZZ_tree_t));

        ZZ_get_data (group, gram, divisors, data, tree, projections, konst_flag);	
	ZZ_intern (gram, data, tree, inzidenz);
	ZZ_put_data (group, data, tree);
	
	if (!GRAPH){
	   /* free tree and data */
	   n = tree->root;
	   do {
              t = n->next;		
	      ZZ_free_node (data, n);
	   } while ((n = t) != NULL);
	   ZZ_free_data(data);
	   free(data);
	   free(tree);
	}
	else{
	   DATEN = (ZZ_super_TYP *)calloc(1,sizeof(ZZ_super_TYP));
	   if (SUBDIRECT) {
		if (PrI) {
			for (i = 0; i < SUBDIRECT; i++) {
				if (PrI[i]) {
					free_mat (PrI[i]);
				}
			}
			free (PrI);
			PrI = NULL;
		}
		if (SUB_VEC) {
			free (SUB_VEC);
			SUB_VEC = NULL;
		}
	   }
	   DATEN->tree = tree;
	   DATEN->data = data;
	}
	
	for (i = 0; i < group->zentr_no; i++) {
		ZZ_transpose_array(group->zentr[i]->array.SZ,
				   group->zentr[i]->cols);
	}
	for (i = 0; i < group->gen_no; i++) {
		ZZ_transpose_array (group->gen[i]->array.SZ, 
				    group->gen[i]->cols);
	}
	
	
	if (ZCLASS == 2 && GRAPH){
	   SUPER_INFO = DATEN;
	}
	if (ZCLASS == 1 && GRAPH){
	   SUPER_info[super_nr] = DATEN;
	}
}













