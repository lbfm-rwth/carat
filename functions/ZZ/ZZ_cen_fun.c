#include <signal.h>
#include "typedef.h"
#include "getput.h"
#include "tools.h"
#include "bravais.h"
#include "base.h"
#include "matrix.h"
#include "longtools.h"
#include "voronoi.h"
#include "ZZ_P.h"
#include "ZZ_zclass_P.h"
#include "ZZ_irr_const_P.h"
#include "ZZ_lll_P.h"
#include "ZZ_cen_fun_P.h"
#include "datei.h"

extern int IDEM_NO;
int *SUB_VEC;
matrix_TYP **PrI;
QtoZ_konst_TYP *KONSTITUENTEN = NULL;

/* eingefuegt von Oliver am 5.2.99 */
int OANZ = 0;
matrix_TYP **OMAT;
int OFLAG = FALSE;	
/* ------------------------------- */


  /*============================================================*\
  ||                                                            ||
  || Enthaelt die speziellen Funktionen des Zentrierungs-       ||
  || algorithmus. Zunaechst die Ein- und Ausgaberoutinen        ||
  || und ihre Hilfsfunktionen.                                  ||
  ||                                                            ||
  \*============================================================*/

/*-------------------------------------------------------------------*\
| takes the matrices in MAT as a generating set for a matrix-ring     |
| over p and evaluates any linearcombination to get the whole ring    |
\*-------------------------------------------------------------------*/

static matrix_TYP **gen_sys (int *num, int l, matrix_TYP ** MAT);

static matrix_TYP **gen_sys (num, l, MAT)
     int *num;
     int l;
     matrix_TYP **MAT;
{
	int i, j, *s, n, anz, w;
	int **Z, **Q1, **Q2;
	int row;
	matrix_TYP **ring;
	int a = 0;

#if 0
	if (l == 1) { /* there`s only the identical Endomorphism and all
			 GF(p) multiples */
		ring = MAT;
		*num = act_prime - 1;
		return (ring);
	}
#endif

	a = intpow (act_prime, l);
	anz = a - 1;
	row = MAT[0]->rows;
	
	ring = (matrix_TYP **) malloc (anz * sizeof (matrix_TYP *));
	*num = anz;
	s = (int *)calloc (l, sizeof (int));
	n = 0;    /* number of the actuell matrix */
	anz = 0;  /* highest coefficient in s     */
	s[0] = 1;
	ring[0] = MAT[0];
	while (anz < l) {
		a = 0;	  /* Laufindex ueber s */
		while ((a <= anz) && ((s[a] = S (s[a], 1)) == 0)) {
			a++;
		}
		if (a <= anz) {
			n++;
			ring[n] = init_mat(row, row, "p");
			Z = Q1 = ring[n]->array.SZ;
			for (a = 0; a <= anz; a++) {
				if ((w = s[a]) == 0) {
					continue;
				}
				Q2 = MAT[a]->array.SZ;
				for (i = 0; i < row; i++) {
					for (j = 0; j < row; j++) {
						Z[i][j] = S(Q1[i][j],
							    P(Q2[i][j], w));
					}
				}
			}
			Check_mat(ring[n]);
		} else {
			if (++anz < l) {
				s[anz] = 1;
				s[anz - 1] = 0;
				n++;
				ring[n] = MAT[anz];
			}
		}
	}
	free (MAT);
	free (s);
	return (ring);
}

void ZZ_free_node (data, n)
     ZZ_data_t *data;
     ZZ_node_t *n;
{
	int i,j;
	ZZ_couple_t *m, *tm;

	if (n->group != NULL){
		free_bravais(n->group);
	}
	if (n->stab_chain != NULL){
		for (i=0;i<n->col_group->dim;i++){
			free_bahn(n->stab_chain[i]);
			free(n->stab_chain[i]);
		}
		free(n->stab_chain);
	}
	if (n->col_group != NULL){
		free_bravais(n->col_group);
	}
	if (n->brav != NULL){
		free_bravais(n->brav);
	}
	if (n->perfect != NULL && n->perfect_no > 0){
		for (i=0;i<n->perfect_no;i++){
			clear_voronoi(n->perfect[i]);
			free(n->perfect[i]);
		}
		free(n->perfect);
		n->perfect = NULL;
		n->perfect_no = 0;
	}
	if (n->N_no_orbits > 0 && n->N_orbits != NULL){
		for (i=0;i<n->N_no_orbits;i++){
			for (j=0;j<n->N_lengths[i];j++)
		 		free_mat(n->N_orbits[i][j]);
			free(n->N_orbits[i]);
		}
		free(n->N_lengths);
		free(n->N_orbits);
	}
	if (n->U != NULL) {
		free_mat (n->U);
	}
	if (n->Q != NULL) {
		free_mat (n->Q);
	}
	if (n->U_inv != NULL) {
		free_mat (n->U_inv);
	}
	if (n->el_div != NULL) {
		free_mat (n->el_div);
	}
	free (n->path);
	for (i = 0; i < data->p_consts.k; i++) {
		free (n->k_vec[i]);
	}
	free (n->k_vec);
	
	if ((m = n->parent) != NULL) {
		do {
			tm = m->elder;
			free (m);
			m = tm;
		} while (m != NULL);
	}
	if ((m = n->child) != NULL) {
		do {
			tm = m->elder;
			free (m);
			m = tm;
		} while (m != NULL);
	}
	free (n);
}

/*----------------------------------------------------------------*\
| calculate the Basis of the endomorphisms for any constituent     |
| and generate the endomorphism-ring. The matrix set TMAT is       |
| copied into data->Endo[i][j] by gen_sys, so there mustn`t be a    |
| ZZ_free_mat at TMAT.                                                |
\*----------------------------------------------------------------*/
void ZZ_make_endo (data)
     ZZ_data_t *data;
{
	int i, j, l;
	matrix_TYP **TMAT;

	for (i = 0; i < data->p_consts.k; i++) {
		init_prime (data->p_consts.p[i]);
		for (j = 0; j < data->p_consts.s[i]; j++) {
			l = data->r;
			TMAT = p_solve(&l,
				       data->p_consts.Delta[i][j],
				       data->p_consts.Delta[i][j], 2);
			data->EnCo[i][j].low = l;
			data->Endo[i][j] = gen_sys(&data->EnCo[i][j].hi,
						   l, TMAT);
		}
	}
}



/*  tests if two constituents are isomorphic. Returns 1 if they are,
 *  0 otherwise.
 *  
 */

static int test_two_constituents(data, i, j, k)
     ZZ_data_t *data;
     int i, j, k;
{
	int ll, l, s;
	matrix_TYP **temp;

#define P_CONST data->p_consts.Delta
	if (data->n[i][j] == data->n[i][k]) {
		for (ll = 0; ll < data->r; ll++) {
			s = 0;
			for (l = 0; l < data->n[i][j]; l++) {
				s = S(S(s, P_CONST[i][j][ll]->array.SZ[l][l]),
				      -P_CONST[i][k][ll]->array.SZ[l][l]);
			}
			if (s != 0) {
				break;
			}
		}
		if (s == 0) {
			s = data->r;
			temp = p_solve(&s, P_CONST[i][j], P_CONST[i][k], 2);
			if (s > 0) {
#if DEBUG
				printf("p: %d, Konst. %d/%d\n", i, j ,k);
#endif
				for (ll = 0; ll < s; ll++) {
					free_mat(temp[ll]);
				}
				free(temp);
				data->p_consts.s[i]--;
				for (ll = 0; ll < data->r; ll++) {
					free_mat(P_CONST[i][k][ll]);
				}
				free(P_CONST[i][k]);
				P_CONST[i][k]= P_CONST[i][data->p_consts.s[i]];
				data->n[i][k]= data->n[i][data->p_consts.s[i]];
				return 1;
			}
		}
	}
	return 0;
#undef P_CONST
}

/*----------------------------------------------------------------*\
| Test if the set of constituents is redundant, i.e two of them    |
| are isomorphic.                                                  |
\*----------------------------------------------------------------*/
void ZZ_test_konst (data)
     ZZ_data_t *data;
{
	int i, j, k;
	
	for (i = 0; i < data->p_consts.k; i++) {
		init_prime (data->p_consts.p[i]);
		for (j = 0; j < data->p_consts.s[i] - 1; j++) {
			for (k = j + 1; k < data->p_consts.s[i]; k++) {
				if (test_two_constituents(data, i, j, k)) {
					k--;
				}
			}
		}			/* Konstituenten */
	}				/* Primzahlen */
}

/*
 *  Diese Funktion fuellt die Datenstrukturen "data" und "tree"
 *  insbesondere werden die irreduziblen Konstituenten fuer die
 *  Primteiler der Gruppenordnung errechnet
 *  q2z ruft diese Fufnktion ueber ZZ mehrmals auf; einige Daten brauchen nicht mehrmals
 *  ausgerechnet werden; das kann noch verbessert werden!
 */
void ZZ_get_data (group, gram, divisors, data, tree, projections, konst_flag)
     bravais_TYP *group;
     matrix_TYP *gram;
     int *divisors;
     ZZ_data_t *data;
     ZZ_tree_t *tree;
     int *projections;
     int konst_flag;
{
	int i, j, k;
	matrix_TYP **help2;
	QtoZ_konst_TYP *data_neu;

	
	if (konst_flag == -3 && GRAPH){
	   for (i = 0; i < KONSTITUENTEN->k; i++){
	      for (j = 0; j < KONSTITUENTEN->s[i]; j++){
	         for (k = 0; k < KONSTITUENTEN->r; k++){
	            free_mat(KONSTITUENTEN->Delta[i][j][k]);
	         }
	         free(KONSTITUENTEN->Delta[i][j]);
	      }
	      free(KONSTITUENTEN->Delta[i]);
	   }
	   free(KONSTITUENTEN->Delta);
	   free(KONSTITUENTEN->s);
	   free(KONSTITUENTEN);	
	   return;
	}
	
	tree->root = tree->last = (ZZ_node_t *) malloc(sizeof(ZZ_node_t));
	data->N = group->dim;

	if (ZCLASS > 0){
		group->gen_no -= IDEM_NO;
		tree->root->col_group = tr_bravais(group,1,FALSE);
		group->gen_no += IDEM_NO;
		tree->root->group = tr_bravais(tree->root->col_group,1,FALSE);
		tree->root->brav = NULL;
		tree->root->stab_chain = NULL;
		tree->root->N_orbits = NULL;
		tree->root->N_lengths = NULL;
		tree->root->N_no_orbits = 0;
		tree->root->perfect = NULL;
		tree->root->perfect_no = 0;
		tree->root->Q = init_mat(data->N, data->N, "1");
        } else {
		tree->root->col_group = NULL;
		tree->root->group = NULL;
		tree->root->brav = NULL;
		tree->root->stab_chain = NULL;
		tree->root->N_orbits = NULL;
		tree->root->N_lengths = NULL;
		tree->root->N_no_orbits = 0;
		tree->root->perfect = NULL;
		tree->root->perfect_no = 0;
		tree->root->Q = NULL;
	}

	tree->root->U = init_mat (data->N, data->N, "");
	tree->root->U_inv = init_mat (data->N, data->N, "");
	/* changed by oliver (6.11.00) from
	tree->root->el_div = init_mat (1, data->N, "");
	to: */
	tree->root->el_div = init_mat (data->N, data->N, "");	
	for (i = 0; i < data->N; i++) {
		tree->root->U->array.SZ[i][i] =
			tree->root->U_inv->array.SZ[i][i]  =
			/* changed by oliver (6.11.00) from
			tree->root->el_div->array.SZ[0][i] = 1;
			to: */
			tree->root->el_div->array.SZ[i][i] = 1;			
	}
	tree->root->number =
		tree->node_count  =
		tree->root->level =
		tree->root->index = 1;
	tree->root->path = (ZZ_prod_t *) calloc(1, sizeof(ZZ_prod_t));
	tree->root->parent = NULL;
	tree->root->next = NULL;
	tree->root->child = NULL;

	/*
	 * Read the number of generators and the generators.
	 */

	data->r = group->gen_no;
	data->DELTA   = (matrix_TYP **)malloc(data->r * sizeof(matrix_TYP *));
	data->DELTA_M = (matrix_TYP **)malloc(data->r * sizeof(matrix_TYP *));
	for (i = 0; i < data->r; i++) {
		data->DELTA[i]   = group->gen[i];
		data->DELTA_M[i] = copy_mat(data->DELTA[i]);
	}

	/*
	 * Read the number of primes
	 */
	for (i = 0, data->p_consts.k = 0; i < 100; i++) {
		if (divisors[i] != 0) {
			data->p_consts.k++;
		}
	}
	data->p_consts.p = (int *)malloc(data->p_consts.k * sizeof (int));
	for (i = 0, j = 0; j < data->p_consts.k; i++) {
		if (divisors[i] != 0) {
			data->p_consts.p[j++] = i;
		}
	}
	data->p_consts.s = (int *) malloc(data->p_consts.k * sizeof(int));
	data->n = (int **) malloc(data->p_consts.k * sizeof(int *));
	data->EnCo = (ZZ_prod_t **)malloc(data->p_consts.k *
					  sizeof(ZZ_prod_t *));
	data->p_consts.Delta =
		(matrix_TYP ****)malloc(data->p_consts.k *
					sizeof (matrix_TYP ***));
	data->Endo =
		(matrix_TYP ****)malloc (data->p_consts.k *
					 sizeof (matrix_TYP ***));
	/*
	 * Read the i'th prime and the number of constituents for i
	 */
	
	/* if - else eingefuegt von Oliver: 8.8.00 */
        if (konst_flag == 0 && GRAPH){
	   for (i = 0; i < KONSTITUENTEN->k; i++) {
              data->p_consts.s[i] = KONSTITUENTEN->s[i];
	      data->n[i] = (int *)malloc(data->p_consts.s[i] * sizeof (int));
	      data->p_consts.Delta[i] =
		   (matrix_TYP ***)malloc(data->p_consts.s[i] *
				          sizeof (matrix_TYP **));
	      for (j = 0; j < data->p_consts.s[i]; j++) {
			   data->p_consts.Delta[i][j] =
				   (matrix_TYP **)malloc(data->r *
						         sizeof (matrix_TYP *));
			   for (k = 0; k < KONSTITUENTEN->r; k++) {
				   data->p_consts.Delta[i][j][k] =
				                copy_mat(KONSTITUENTEN->Delta[i][j][k]);
				   data->p_consts.Delta[i][j][k]->prime =
					   data->p_consts.p[i];
			   }
			   data->n[i][j] = data->p_consts.Delta[i][j][0]->rows;		
	      }
           }
        }
        else{
           /* alte Version */
	   for (i = 0; i < data->p_consts.k; i++) {
	   	   help2 = ZZ_irr_const (data->DELTA, data->r,
				         data->p_consts.p[i],
				         &data->p_consts.s[i]);

		   data->n[i] = (int *)malloc(data->p_consts.s[i] * sizeof (int));
		   data->p_consts.Delta[i] =
			   (matrix_TYP ***)malloc(data->p_consts.s[i] *
					          sizeof (matrix_TYP **));
		   for (j = 0; j < data->p_consts.s[i]; j++) {
			   data->p_consts.Delta[i][j] =
				   (matrix_TYP **)malloc(data->r *
						         sizeof (matrix_TYP *));
			   for (k = 0; k < data->r; k++) {
				   data->p_consts.Delta[i][j][k] =
					   help2[j * data->r + k];
				   data->p_consts.Delta[i][j][k]->prime =
					   data->p_consts.p[i];
			   }
			   data->n[i][j] = data->p_consts.Delta[i][j][0]->rows;
		   }
		   free (help2);
	   }
	   ZZ_test_konst(data);
	   ZZ_test_konst(data);
	   ZZ_test_konst(data);
	}
	
	/*
	 * initialize Endomorphisms
	 */
	for (i = 0; i < data->p_consts.k; i++) {
		data->EnCo[i] =
			(ZZ_prod_t *)malloc(data->p_consts.s[i] *
					    sizeof (ZZ_prod_t));
		data->Endo[i] =
			(matrix_TYP ***)malloc(data->p_consts.s[i] *
					       sizeof (matrix_TYP **));
	}

	data->epi_base = NULL;
	data->epi = init_mat(data->N, data->N, "");
	tree->root->k_vec = (int **)malloc(data->p_consts.k * sizeof(int *));
	data->VK = (int **) malloc (data->p_consts.k * sizeof (int *));
	for (i = 0; i < data->p_consts.k; i++) {
		tree->root->k_vec[i] =
			(int *)calloc(data->p_consts.s[i], sizeof (int));
		data->VK[i] = (int *)calloc(data->p_consts.s[i] + 1,
					    sizeof(int));
		data->VK[i]++;
	}
	ZZ_make_endo (data);
	if (SUBDIRECT) {
		SUBDIRECT = projections[0];
		SUB_VEC = (int *) malloc (SUBDIRECT * sizeof (int));
#ifdef DEBUG
		fprintf (stderr, "SUBDIRECT = %d\n", SUBDIRECT);
#endif
		for (i = 0; i < SUBDIRECT; i++) {
			SUB_VEC[i] = projections[i+1];
#ifdef DEBUG
			fprintf(stderr, "Dimension[%d] = %d\n", i, SUB_VEC[i]);
#endif
		}
		PrI = (matrix_TYP **) malloc(SUBDIRECT * sizeof(matrix_TYP *));
		k = 0;
		for (i = 0; i < SUBDIRECT; i++) {
			PrI[i] = init_mat (data->N, SUB_VEC[i], "");
			for (j = 0; j < SUB_VEC[i]; j++) {
				PrI[i]->array.SZ[k + j][j] = 1;
#ifdef DEBUG
				fput_mat (stderr, PrI[i], "PrI[i]", 0);
#endif
			}
			k += SUB_VEC[i];
		}
	}
	
	
	if (konst_flag == 1 && GRAPH){
	   KONSTITUENTEN = (QtoZ_konst_TYP *)calloc(1, sizeof(QtoZ_konst_TYP));
	   KONSTITUENTEN->k = data->p_consts.k;
	   KONSTITUENTEN->s = (int *)calloc(data->p_consts.k, sizeof(int));
	   KONSTITUENTEN->r = (data->r - IDEM_NO);
	   KONSTITUENTEN->Delta = (matrix_TYP ****)malloc(KONSTITUENTEN->k *
					sizeof (matrix_TYP ***));

	   for (i = 0; i < data->p_consts.k; i++) {
              KONSTITUENTEN->s[i] = data->p_consts.s[i];
  	      KONSTITUENTEN->Delta[i] = (matrix_TYP ***)malloc(KONSTITUENTEN->s[i] *
					          sizeof (matrix_TYP **));
		   for (j = 0; j < KONSTITUENTEN->s[i]; j++) {
			   KONSTITUENTEN->Delta[i][j] = (matrix_TYP **)malloc(KONSTITUENTEN->r *
						         sizeof (matrix_TYP *));
			   for (k = 0; k < KONSTITUENTEN->r; k++) {
				   KONSTITUENTEN->Delta[i][j][k] =
					   copy_mat(data->p_consts.Delta[i][j][k]);
			   }
		   }
	   }
	}
}

/*****************************************************************************/
/*                                                                           */
/* input: global data-structure in data and the tree of Zentrierungen        */
/* output: the structure-component group->zentr points to an array of          */
/*         of matrices, each containing a set of generators of a             */
/*         Zentrierung as Z-modul. group->zentr_no contains the number         */
/*         of the computed Zentrierungen                                     */
/*                                                                           */
/*****************************************************************************/
void ZZ_put_data (group, data, tree)
     bravais_TYP *group;
     ZZ_data_t *data;
     ZZ_tree_t *tree;
{
	ZZ_node_t *n, *t;
	int i;
	ZZ_couple_t *m;
	
	group->zentr = (matrix_TYP **)malloc(tree->node_count *
					     sizeof (matrix_TYP *));
	group->zentr_no = tree->node_count;
	n = tree->root;
	i = 0;
	do {
		group->zentr[i] = copy_mat(n->U);
		/* n->U = NULL; oliver: 16.8.00 */
 		if (SHORTLIST) {
			if (SUBDIRECT) {
				printf ("\t L%d: Number of Sublattices: %d ",
					n->number, n->anz_tg);
			} else {
				printf ("\t L%d: ", n->number);
			}
			if ((m = n->child) != NULL) {
				do {
					printf (" L%d ", m->he->number);
				}
				while ((m = m->elder) != NULL);
			}
			printf ("\n");
			fflush (stdout);
		}
		i++;
		t = n->next;		
		/* oliver: 11.8.00
		ZZ_free_node (data, n); */
	} while ((n = t) != NULL);
}

void ZZ_fput_data (data, tree, ABBRUCH)
     ZZ_data_t *data;
     ZZ_tree_t *tree;
     int ABBRUCH;
{
	ZZ_node_t *n, *t;
	ZZ_couple_t *m;
	int old_lev = 0, lev_cnt = 0;

	n = tree->root;
	if (TEMPORAER) {
		fprintf (ZZ_temp,"\n==========================================================================\n");
	}
	do {
		if (old_lev != n->level) {
			if (TEMPORAER && (old_lev)) {
				fprintf (ZZ_temp, "\n======= %d Lattices on level %d ====================================\n", lev_cnt, old_lev);
			}
			old_lev = n->level;
			lev_cnt = 0;
		}
		if (TEMPORAER && !(ABBRUCH & 1)) {
			fprintf (ZZ_temp,
				 "\n\t      L%3.d: %13s\tConstituent\n",
				 n->number, "Lattice");
			if (SUBDIRECT) {
				fprintf (ZZ_temp,
					 "\n Number of sublattices: %d ",
					 n->anz_tg);
			}
			if ((m = n->parent) != NULL) {
				fprintf (ZZ_temp, "     parents  :\n");
				do {
					fprintf(ZZ_temp, 
						"%25s%-3d\t(p=%d,Nr.=%d)\n",
						"L", m->he->number,
						data->p_consts.p[m->she.i],
						m->she.j + 1);
				} while ((m = m->elder) != NULL);
			}
			if ((m = n->child) != NULL) {
				fprintf (ZZ_temp, "     children :\n");
				do {
					fprintf (ZZ_temp,
						 "%25s%-3d\t(p=%d,Nr.=%d)\n",
						 "L",
						 m->he->number,
						 data->p_consts.p[m->she.i],
						 m->she.j + 1);
				} while ((m = m->elder) != NULL);
			}
		}
		lev_cnt++;
		t = n->next;
	} while ((n = t) != NULL);
	if (TEMPORAER) {
		fprintf (ZZ_temp,
			 "\n======= %d Lattices on level %d ====================================\n", lev_cnt, old_lev);
	}
	if (ABBRUCH & 2) {
		if (TEMPORAER) {
			fprintf (ZZ_temp, " \t ============================================================\n \t Algorithmus nach Berechnung von %d Zentrierungen abgebrochen \n \t ============================================================\n", NUMBER);
		}
	}
	if (ABBRUCH & 4) {
		if (TEMPORAER) {
			fprintf (ZZ_temp, " \t ============================================================\n \t Algorithmus nach Berechnung von %d Levels abgebrochen \n \t ============================================================\n", LEVEL);
		}
	}
}

void ZZ_free_data (data)
     ZZ_data_t *data;
{
	int i, j, k;
	
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
	
	for (i = 0; i < data->p_consts.k; i++) {
		for (j = 0; j < data->p_consts.s[i]; j++) {
			for (k = 0; k < data->r; k++) {
				free_mat (data->p_consts.Delta[i][j][k]);
			}
			for (k = 0; k < data->EnCo[i][j].hi; k++) {
				free_mat (data->Endo[i][j][k]);
			}
			free (data->p_consts.Delta[i][j]);
			free (data->Endo[i][j]);
		}
		free (data->p_consts.Delta[i]);
		free (data->Endo[i]);
		free (data->EnCo[i]);
		data->VK[i]--;
		free (data->VK[i]);
		free (data->n[i]);
	}
	
	free (data->p_consts.Delta);
	free (data->Endo);
	free (data->EnCo);
	free (data->VK);
	free (data->n);
	
	for (i = 0; i < data->r; i++) {
		free_mat (data->DELTA_M[i]);
	}
	free (data->DELTA);
	free (data->DELTA_M);
	
	free_mat (data->epi);
	free (data->p_consts.s);
	free (data->p_consts.p);
	if (ZZ_temp) {
		fclose (ZZ_temp);
		ZZ_temp = NULL;
	}
}

/*------------------------------------------------------------*\
| Variation of p_solve, it differs in the way of output        |
\*------------------------------------------------------------*/
matrix_TYP *ZZ_p_vsolve (anz, L_mat, R_mat)
     int *anz;
     matrix_TYP **L_mat, **R_mat;

{
	int **M, **N, *v, f, help;
	int *M_j, *M_act, numax;
	matrix_TYP *E;
	boolean L_null, R_null;
	int *ss, *nuin, p, i, j, jj, k, kk, l, cR, cL, cM, rR, rL, rM, rN;
	int rg = 0, act;

	
	p = act_prime;
	L_null = (L_mat == NULL);
	R_null = (R_mat == NULL);
	
	cL = L_null ? 1 : L_mat[0]->cols;
	rL = L_null ? 1 : L_mat[0]->rows;
	rR = R_null ? 1 : R_mat[0]->rows;
	cR = R_null ? 1 : R_mat[0]->cols;
	
	cM = cL * rR;
	rM = *anz * rL * cR;
	rN = rL * cR;
	rg = *anz * cR;
	N = (int **)malloc(rN * sizeof(int *));
	M = (int **)malloc(rM * sizeof(int *));

	/* changed tilman 6/2/97 from
	 * nuin = (int *)malloc(cM * sizeof(int));
	 * to
	 */
	nuin = (int *)malloc((cM+1) * sizeof(int));

	ss = (int *) malloc(cM * sizeof(int));

	for (i = 0; i < *anz; i++) {
		if (!L_null) {
			if ((L_mat[i]->cols != cL) || (L_mat[i]->rows != rL)) {
				fprintf (stderr, "Error in input\n");
				exit (3);
			}
		}
		if (!R_null) {
			if ((R_mat[i]->rows != rR) || (R_mat[i]->cols != cR)) {
				fprintf (stderr, "Error in input\n");
				exit (3);
			}
		}
		
		for (j = 0; j < rN; j++) {
			N[j] = (int *) calloc (cM, sizeof (int));
		}
		
		if (!L_null) {
			for (j = 0, jj = 0; j < rL; j++, jj += rR) {
				for (k = 0, kk = 0; k < cL; k++, jj -= rR) {
					if ((help=L_mat[i]->array.SZ[j][k] % p)
					    < 0) {
						help += p;
					}
					if (help) {
						for (l = 0;
						     l < rR; l++, jj++, kk++) {
							N[jj][kk] = help;
						}
					} else {
						jj += rR;
						kk += rR;
					}
				}
			}
		}
		if (!R_null) {
			for (j = 0; j < cR; j++) {
				for (k = 0; k < rR; k++) {
					if ((help=R_mat[i]->array.SZ[k][j] % p)
					    < 0) {
						help += p;
					}
					if (help) {
						for (l = 0, jj = j, kk = k;
						     l < cL;
						     l++, jj += cR, kk += rR) {
							N[jj][kk]= S(N[jj][kk],
								     -help);
						}
					}
				}
			}
		}
		
		/* permutierte Eingabe in Liste */
		for (j = 0; j < rL; j++) {
			for (k = 0; k < cR; k++) {
				M[j*rg + i*cR + k] = N[(rL - 1 - j) * cR + k];
			}
		}
	}
	rg = 0;
	/*
	 * Gauss Algorithm
	 */
	for (i = 0; i < cM; i++) {
		ss[i] = -1;
		/*
		 * Find the row with non-zero entry in i-th column
		 */
		for (j = rg; (j < rM) && (M[j][i] == 0); j++);
		if (j == rM) {
			continue;
		}
		act = j;
		/*
		 * Normalize act-th row and mark the non-zero entries
		 */
		nuin[0] = 1;
		f = P (1, -M[j][i]);
		for (k = i; k < cM; k++) {
			if ((help = M[act][k])) {
				M[act][k] = P (help, f);
				nuin[nuin[0]++] = k;
			}
		}
		/*
		 * Swap act-th row and rg-th row
		 */
		v = M[rg];
		M[rg] = M[act];
		M[act] = v;
		ss[rg] = i;
		act = rg++;
		/*
		 *  Clear i-th column  downwards
		 */
		M_act = M[act];
		for (j = act + 1; j < rM; j++) {
			if ((f = S (0, -M[j][i])) != 0) {
				M_j = M[j];
				M_j[i] = 0;
				numax = nuin[0];
				if (f == 1) {
					for (k = 2, kk = nuin[2];
					     k < numax; kk = nuin[++k]) {
						M_j[kk] = S(M_j[kk],M_act[kk]);
					}
				} else {
					for (k = 2, kk = nuin[2];
					     k < numax; kk = nuin[++k]) {
						M_j[kk] = S(M_j[kk], 
							    P(M_act[kk], f));
					}
				}
			}
		}
	}
	if ((*anz = cM - rg) != 0) {
		/*
		 * Matrix has not full rank: clear it upwards
		 */
		for (i = rg - 1; i > 0; i--) {
			nuin[0] = 2;
			nuin[1] = ss[i];
			M_act = M[i];
			for (j = ss[i] + 1; j < cM; j++) {
				if (M_act[j]) {
					nuin[nuin[0]++] = j;
				}
			}
			if (nuin[0] == 2) {
				j = ss[i];
				for (k = i - 1; k > 0; k--) {
					M[k][j] = 0;
				}
			} else {
				for (j = i - 1; j >= 0; j--) {
					if ((f = S (0, -M[j][ss[i]])) != 0) {
						M_j = M[j];
						M_j[ss[i]] = 0;
						numax = nuin[0];
						if (f == 1) {
							for (k= 2,kk = nuin[2];
							     k < numax;
							     kk = nuin[++k]) {
					M_j[kk] = S(M_j[kk], M_act[kk]);
							}
						} else {
							for (k= 2, kk= nuin[2];
							     k < numax;
							     kk = nuin[++k]) {
				M_j[kk] = S (M_j[kk], P (M_act[kk], f));
							}
						}
					}
				}
			}
		}
		E = init_mat (rR, rR, "");
		for (i = 0, k = 0, l = 0; i < cM; i++) {
			if (i == ss[k]) {
				E->array.SZ[l + (i / rR)][i % rR] = p;
				l++;
				k++;
				continue;
			}
			for (j = 0; j < k; j++) {
				if (M[j][i] != 0) {
					E->array.SZ
						[l + (ss[j] / rR)]
						[ss[j] % rR] = p - M[j][i];
				}
			}
			E->array.SZ[l + (i / rR)][i % rR] = 1;
			l++;
		}
	}
	for (i = 0; i < rM; i++) {
		free (M[i]);
	}
	free (M);
	free (N);
	free (ss);
	free (nuin);
	
	return (E);
}

/*---------------------------------------------------------------*\
| Take the kk-th epimorphism from DELTA_M onto p_consts.Delta[ii][jj],    |
| determine its kernel, extend to a generating set of the new    |
| centering and calculate the common factors and index of that    |
| centering. Return the possibly new node for checking, wether it |
| occured already earlier in the algorithm or not.          |
\*---------------------------------------------------------------*/
ZZ_node_t *ZZ_center (data, father, ii, jj)
     ZZ_data_t *data;
     ZZ_node_t *father;
     int ii, jj;
{
	int flag, sg, i, j, d;
	matrix_TYP *ker;
	ZZ_node_t *n;
	
	n = (ZZ_node_t *) malloc (sizeof (ZZ_node_t));
	n->col_group = NULL;
	n->group = NULL;
	n->brav = NULL;
	n->N_orbits = NULL;
	n->N_lengths = NULL;
	n->N_no_orbits = 0;
	n->stab_chain = NULL;
	n->perfect = NULL;
	n->perfect_no = 0;
	n->number = -1;
	n->index = -1;
	n->level = father->level + 1;
	n->U = n->Q = n->U_inv = n->el_div = NULL;
	n->next = NULL;
	n->parent = n->child = NULL;
	n->k_vec = (int **) malloc (data->p_consts.k * sizeof (int *));
	for (i = 0; i < data->p_consts.k; i++) {
		n->k_vec[i] =
			(int *)malloc(data->p_consts.s[i] * sizeof (int));
		memcpy (n->k_vec[i], father->k_vec[i],
			data->p_consts.s[i] * sizeof (int));
	}
	n->k_vec[ii][jj]++;
	d = father->level;
	n->path = (ZZ_prod_t *) malloc ((d + 1) * sizeof (ZZ_prod_t));
	memcpy (n->path, father->path, d * sizeof (ZZ_prod_t));
	n->path[d].hi = ii;
	n->path[d].low = jj;
	
	/*------------------------------------------------------------*\
	  | Determine the kernel of the epimorphism            |
	  \*------------------------------------------------------------*/
	d = 1;
	ker = ZZ_p_vsolve (&d, NULL, &data->epi);
	if (d != 0) {
		/*
		 * Express that kernel in the basis of the original lattice
		 */
		n->U = mat_mul (ker, father->U);
		free_mat (ker);
		/*
		 * Make n->U a generating set of the new ZZ_centering by
		 * appending to n->U  p times the basis of father.
		 */
	} else {
		n->U = copy_mat (father->U);
		iscal_mul (n->U, act_prime);
	}
	/*
	 * Calculate the invariant factors of that generating set
	 */
	Check_mat (n->U);
	
	if (U_option) {
		n->el_div = long_elt_mat(NULL,n->U, NULL);
		/*
		 * Calculate index as product of the invariant factors
		 */
		n->index = n->el_div->array.SZ[0][0];
		for (i = 1; i < n->el_div->cols; i++) {
		        /* changed by oliver (6.11.00) from:
			n->index *= n->el_div->array.SZ[0][i];
			to: */
			n->index *= n->el_div->array.SZ[i][i];
		}
	} else {
		n->el_div = init_mat (1, 1, "");
		sg = n->U->array.SZ[0][0];
		for (i = 0; ((sg != 1) && (i < n->U->rows)); i++) {
			for (j = 0; ((sg != 1) && (j < n->U->cols)); j++) {
				sg = GGT (sg, n->U->array.SZ[i][j]);
			}
		}
		n->el_div->array.SZ[0][0] = (sg > 0) ? sg : -sg;
	}
	sg = n->el_div->array.SZ[0][0];
	if (sg != 1) {
		flag = TRUE;
		for (i = 0; flag && (i < data->p_consts.k); i++) {
			if (data->p_consts.p[i] == sg) {
				if (data->VK[i][-1] != 0) {
					for (j = 0;
					     j < data->p_consts.s[i]; j++) {
						n->k_vec[i][j]-=data->VK[i][j];
					}
					n->level -= data->VK[i][-1];
				} else {
					data->VK[i][-1] = n->level - 1;
					memcpy (data->VK[i], n->k_vec[i],
						data->p_consts.s[i] *
						sizeof (int));
					memset(n->k_vec[i], 0,
					       data->p_consts.s[i] *
					       sizeof (int));
					n->level = 1;
				}
				flag = FALSE;
			}
		}
	}
	return (n);
}

/*---------------------------------------------------------------*\
| Let E be an epimorphism from DELTA_M onto p_consts.Delta[ii][jj]      |
| Solve the matrix-equation DELTA_M * E = E * p_consts.Delta[ii][jj]     |
| over F(p) and return a basis of the space of all solutions    |
\*---------------------------------------------------------------*/
int ZZ_epimorphs (data, ii, jj)
     ZZ_data_t *data;
     int ii, jj;
{
	int d, i, j, k, x, e;
	int **X, **M;
	int *bas, found, row, col, sae;
	int actlist, vor;
	int a = 0;
	matrix_TYP *TMP, *AMP, **liste;
	
	init_prime (data->p_consts.p[ii]);
	
	d = data->r;
	data->epi_base = p_solve(&d, data->DELTA_M,
				 data->p_consts.Delta[ii][jj], 2);
	if (d > 1) {
		data->epi->cols = data->epi_base[0]->cols;
		if (data->EnCo[ii][jj].low != 1) {
			a = intpow (act_prime, d);
			if (d == data->EnCo[ii][jj].low) {
				/* All Epimorphisms differs only by an
				 * endomorphism 
				 */
				for (i = 1; i < d; i++)	{
					free_mat (data->epi_base[i]);
				}
				data->epi_base =
				  (matrix_TYP **)realloc(data->epi_base,
							 sizeof(matrix_TYP *));
				d = 1;
				return (d);
			}
			sae = d / data->EnCo[ii][jj].low;
			liste = (matrix_TYP **) malloc(a*sizeof(matrix_TYP *));
			actlist = 0;
			bas = (int *) calloc (d, sizeof (int));
			row = data->epi_base[0]->rows;
			col = data->epi_base[0]->cols;
			found = 0;
			for (x = 0; x < d - 1; x++) {
				if (bas[x] == 1) { 
					/* epim already found previously */
					continue;
				}
				if ((++found) == sae) {
					/* got a complete set of epi`s */
					for (++x; x < d; x++) {
						bas[x] = 1;
					}
					break;
				}
				/* Calculate the produts of the
				 * Epimorphism X with all
				 * Endomorphisms and compare to the
				 * other Epi`s
				 */
				vor = actlist;
				for (e = 0; e < data->EnCo[ii][jj].hi; e++) {
					TMP = mat_mul (data->epi_base[x],
						       data->Endo[ii][jj][e]);
					liste[actlist++] = TMP;
					M = TMP->array.SZ;
					for (i = x + 1; i < d; i++) {
						X= data->epi_base[i]->array.SZ;
						j = 0;
						while (
		       (j < row) && ((memcmp (M[j], X[j], 4 * col)) == 0)) {
							j++;
						}
						if (j == row) {
							bas[i] = 1;
							break;
						}
					}
					for (k = 0; k < vor; k++) {
						AMP = pmat_add (TMP, liste[k],
								1, 1);
						liste[actlist++] = AMP;
						for (i = x + 1; i < d; i++) {
					X = data->epi_base[i]->array.SZ;
					j = 0;
					while (
		       (j < row) && ((memcmp (M[j], X[j], 4 * col)) == 0)) {
						j++;
					}
					if (j == row) {
						bas[i] = 1;
						break;
					}
						}
					}
				}
			}
			/*
			 * prepare the new epi_base for output
			 */
			for (i = 0; i < actlist; i++) {
				free_mat (liste[i]);
			}
			free (liste);
			found = 1;
			for (i = 1; i < d; i++) {
				if (bas[i] == 1) {
					free_mat (data->epi_base[i]);
				} else {
					data->epi_base[found++] =
						data->epi_base[i];
				}
			}
			d = found;
			data->epi_base =
				(matrix_TYP **)realloc(data->epi_base,
						       d*sizeof(matrix_TYP *));
		}
	} else {
		if (d != 0) {
			data->epi->cols = data->epi_base[0]->cols;
		} else {
			data->epi->cols = 0;
		}
	}
	return (d);
}

boolean ZZ_successor (data, act)
     ZZ_data_t *data;
     ZZ_node_t **act;

{
	int i;
	
	if ((*act = (*act)->next) == NULL) {
		return (FALSE);
	}
	for (i = 0; i < data->r; i++) {
		free_mat (data->DELTA_M[i]);
		data->DELTA_M[i] = mat_kon((*act)->U, data->DELTA[i],
					   (*act)->U_inv);
	}
	return (TRUE);
}



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
int ZZ_ins_node (Gram, data, tree, father, new, ii, jj, inzidenz, nr, NEU, flagge, g, nnn)
     matrix_TYP *Gram;
     ZZ_data_t *data;
     ZZ_tree_t *tree;
     ZZ_node_t *father, *new;
     int ii, jj;
     QtoZ_TYP *inzidenz;
     int *nr;
     int *NEU;
     int *flagge;
     int *g;
     ZZ_node_t **nnn;
{
	ZZ_node_t *n;
	ZZ_couple_t *c;
	matrix_TYP *Tmp1, *Tmp2, *U_vor, *tmp;
	matrix_TYP *Mat, *Trf, *el, *GMat;
	int **U, **CC, sum, f;
	int i, j, k, nl, flag, ff, gg;
	int ABBRUCH;

	/* inserted tilman 6/2/97 */
	char filename[32];

	sum =
		f =
		g[0] = 0;
	U = new->U->array.SZ;
	ABBRUCH = FALSE;

	/*
	 * Determine gcd of common factors of new->U
	 */
	
	if (new->el_div->array.SZ[0][0] == 1) {
		n = father;
	} else {
		n = tree->root;
	}
	while ((n != NULL) && (n->level < new->level)) {
		n = n->next;
	}

	if (n != NULL) {
		/* We found a node on the same level, now look for
                   all on this level for an identical lattice */
		g[0] = new->el_div->array.SZ[0][0];
		nl = new->level;
		do {
			f = n->U_inv->kgv;
			f *= g[0];
			flag = TRUE;
			for (i = 0; (i < data->p_consts.k) && flag; i++) {
				flag = memcmp (n->k_vec[i],
					       new->k_vec[i],
					       data->p_consts.s[i] *
					       sizeof (int)) == 0;
			}
			if (flag) {
				CC = n->U_inv->array.SZ;
				/*
				 * Index-condition satisfied
				 */

				/* Multiply new->U and n->U_inv and check
				 * divisibility condition for each entry
				 */
				for (i = 0; i < new->U->rows; i++) {
					for (j = 0; j < new->U->cols; j++) {
						sum = 0;
						for (k = 0;
						     k < n->U_inv->rows;
						     k++) {
						sum += U[i][k] * CC[k][j];
						}
						if ((sum % f) != 0) {
							goto next;
						}
					}
				}
				c = father->child;
				while (c != NULL) {
					if (c->he == n) {
						ZZ_free_node (data, new);
				                /* next 7 lines by oliver: 9.8.00: for graph for QtoZ */
						if (ZCLASS == 1 && GRAPH){
						   nr[0] = suche_mat(n->U, inzidenz->gitter, inzidenz->anz);
						   if (nr[0] == -1){
						      fprintf(stderr,"ERROR 1 in ZZ_ins_node!\n");
						      exit(2);
						   }
						   flagge[0] += 1;
						}
						return ABBRUCH;
					}
					c = c->elder;
				}
				break;
			}
		next:
			continue;
		} while (((n = n->next) != NULL) && (n->level == new->level));
	}
	
	if ((n == NULL) || (n->level > nl)) {
		/* this is really a new lattice */

		if (SUBDIRECT) {
			for (i = 0; i < SUBDIRECT; i++) {
				Tmp1 = mat_mul (new->U, PrI[i]);
				Tmp2 = long_elt_mat(NULL,Tmp1, NULL);
				free_mat (Tmp1);
				/* changed 30/06/97 tilman form:
				if (Tmp2->array.SZ[0][Tmp2->cols - 1] != 1) {
				to: */
				if (Tmp2->array.SZ[Tmp2->cols - 1]
                                                  [Tmp2->cols - 1] != 1) {
					i = SUBDIRECT + 1;
				}
				free_mat (Tmp2);
			}
			if (i > SUBDIRECT) {
				ZZ_free_node (data, new);
				
				/* next 2 lines by oliver: 9.8.00: for graph for QtoZ */
				if (GRAPH)
				   nr[0] = -1;
				
				return ABBRUCH;
			}
		}

		if (ZCLASS == 1){
		        /* second call of ZZ by q2z */
			if (orbit_under_normalizer(data,tree,father,new,ii,jj,inzidenz,nr,nnn)){
				ZZ_free_node (data, new);
				flagge[0] += 10;
				return ABBRUCH;
			}
                }
		if (ZCLASS == 2){
		        /* first call of ZZ by q2z */
			if (deal_with_ZCLASS(data, tree, father, new)){
				ZZ_free_node (data, new);
				return ABBRUCH;
			}
		}

		/*
		 * New lattice:
		 * Calculate the scalarproducts of the generating system
		 */
		NEU[0] = 1;
		
		if (ZCLASS == 2 && GRAPH){
		   U_vor = mat_inv(new->U);
		}
		
		/* !!!!!!!!!!!!!! new->U is changed in scal_pr !!!!!!!!!!!!!!!!!!! */
		GMat = Mat = scal_pr (new->U, Gram, FALSE);
		if (ZCLASS == 2 && GRAPH){
		   /* save transformation matrix */
		   tmp = mat_mul(new->U, U_vor);
		   new->Q = tr_pose(tmp);
		   free_mat(tmp);
		   free_mat(U_vor);
		}
		
		if (ZCLASS == 1 && GRAPH){
	           /* now calculate the (integral) representation
                      on the new lattice (bare in mind that it is
   	              row invariant. There is a flaw: normalizer and
 	              centralizer won't be correct */
	           ff = tree->root->group->normal_no;
 	           tree->root->group->normal_no = 0;
	           gg = tree->root->group->cen_no;
 	           tree->root->group->cen_no = 0;
	           new->group = konj_bravais(tree->root->group, new->U);
	           tree->root->group->normal_no = ff;
	           tree->root->group->cen_no = gg;

	           /* the second flaw of konj_bravais is the formspace */
	           for (f = 0; f < new->group->form_no; f++)
		      new->group->form[f]->kgv = 1;
	           long_rein_formspace(new->group->form, new->group->form_no, 1);
	           new->col_group = tr_bravais(new->group, 1, FALSE);
	        }
	
		if (LLLREDUCED)	{
			/*
			 * Perform MLLL-Reduction
			 */
			Trf = ZZ_lll (Mat, 0);
			real_mat (Trf, Mat->rows, Trf->cols);
			/* Store the new Gram matrix and the invariant
			 * factors
			 */
			GMat = Mat;
			/*  Apply the LLL-Transformation to the basis
			 */
			Mat = new->U;
			new->U = mat_mul (Trf, Mat);
			free_mat (Mat);
			free_mat(Trf);
		}
		if (G_option) {
			el = long_elt_mat(NULL,GMat, NULL);
		}
		/*
		 * Invert the Basis-Transformation
		 */
		new->U_inv = mat_inv (new->U);
		new->number = ++tree->node_count;
		new->level = father->level + 1;
		tree->last->next = new;
		tree->last = new;
		
		if (new->number > NUMBER) {
			ABBRUCH = 2;
		}
		if (new->level >= LEVEL) {
			ABBRUCH = 4;
		}

		if (!QUIET) {
			fprintf (stderr, "L%d with (%d,%d) yields L%d\n",
				 father->number, ii, jj, tree->node_count);
		}
		/*
		 * write new lattice to temporary file ZZ.tmp
		 */
		if (TEMPORAER && !NURUMF) {
			fprintf (ZZ_temp, "L%d with (%d,%d) yields L%d\n",
				 father->number, ii, jj, tree->node_count);
			
			if (U_option) {
				fput_mat (ZZ_temp, new->el_div,
					  "Elementary divisors    :", 2);
			}
			ZZ_transpose_array(new->U->array.SZ, new->U->cols);
			/* fput_mat (ZZ_temp, new->U,
				  "Change of basis          :", 4); */
			ZZ_transpose_array(new->U->array.SZ, new->U->cols);
			
			if (G_option) {
				fput_mat(ZZ_temp, el,
				 "Elementary divisors of the Gram matrix:", 2);
			}
			fput_mat (ZZ_temp, GMat,
				  "Gram matrix               :", 2);
			fflush (ZZ_temp);
		}
		if (TEMPORAER && NURUMF) {
			ZZ_transpose_array(new->U->array.SZ, new->U->cols);
			fput_mat (ZZ_temp, new->U,
				  "Change of basis:", 4);
			ZZ_transpose_array(new->U->array.SZ, new->U->cols);
			ZZ_transpose_array(new->U_inv->array.SZ, new->U->cols);
			fput_mat (ZZ_temp, new->U_inv,
				  "Inverse of it  :", 4);
			ZZ_transpose_array(new->U_inv->array.SZ, new->U->cols);
			fflush (ZZ_temp);
		}
		
		/* eingefuegt von Oliver am 5.2.99 */
                if (OFLAG){
			OMAT[OANZ] = copy_mat(new->U);
      			OANZ += 1;
                }
		/* ------------------ */
		
		if (G_option) {
			free_mat (el);
		}
		free_mat (GMat);
	} else {
		/*
		 * We got that lattice already
		 */
		if (!QUIET) {
			fprintf (stderr, "L%d with (%d,%d) yields %d.L%d\n",
				 father->number, ii, jj, g[0], n->number);
		}
		if (TEMPORAER && !NURUMF) {
			fprintf (ZZ_temp, "L%d with (%d,%d) yields %d.L%d\n",
				 father->number, ii, jj, g[0], n->number);
		}
		
                /* next 7 lines by oliver: 9.8.00: for graph for QtoZ */
		if (ZCLASS == 1 && GRAPH){		
		   nr[0] = suche_mat(n->U, inzidenz->gitter, inzidenz->anz);
		   if (nr[0] == -1){
		      fprintf(stderr,"ERROR 2 in ZZ_ins_node!\n");
		      exit(3);
		   }
		   flagge[0] += 100;
		}
		ZZ_free_node (data, new);
		new = n;
	}
	
	/*
	 * Tell the son, who's father
	 */
	c = (ZZ_couple_t *) malloc (sizeof (ZZ_couple_t));
	c->he = father;
	c->she.i = ii;
	c->she.j = jj;
	c->factor = g[0];
	c->elder = new->parent;
	new->parent = c;

	/*
	 * Tell father, he's got a new son
	 */
	c = (ZZ_couple_t *) malloc (sizeof (ZZ_couple_t));
	c->he = new;
	c->she.i = ii;
	c->she.j = jj;
	c->factor = g[0];
	c->elder = father->child;
	father->child = c;
	return ABBRUCH;
}




/*}}}  */
/*{{{  ZZ_pick_epi */
void 
ZZ_pick_epi (data, number, ii, jj)
     ZZ_data_t *data;
     int number, ii, jj;
{
  int **E, **hh;
  int q, i, j, k, f;
  matrix_TYP *H;

  init_prime (data->p_consts.p[ii]);
  q = data->EnCo[ii][jj].hi + 1;
  E = data->epi->array.SZ;
  for (i = 0; i < data->epi->rows; i++)
    {
      for (j = 0; j < data->epi->cols; j++)
	{
	  E[i][j] = 0;
	}
    }
  for (i = 0; number != 0; number /= q, i++)
    {
      if ((f = (number % q)) != 0)
	{
	  H = mat_mul (data->epi_base[i], data->Endo[ii][jj][f - 1]);
	  hh = H->array.SZ;
	  for (j = 0; j < data->epi->rows; j++)
	    {
	      for (k = 0; k < data->epi->cols; k++)
		{
		  E[j][k] = S (E[j][k], hh[j][k]);
		}
	    }
	  free_mat (H);
	}
    }
}

/*}}}  */
