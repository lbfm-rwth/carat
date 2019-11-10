#ifndef CARAT_ZZ_CEN_FUN_P_H
#define CARAT_ZZ_CEN_FUN_P_H

#include "typedef.h"
#include "matrix.h"
#include "ZZ_P.h"

extern void ZZ_free_node(ZZ_data_t * data, ZZ_node_t * n);
extern void ZZ_make_endo(ZZ_data_t * data);
extern void ZZ_test_konst(ZZ_data_t * data);
extern void ZZ_get_data(bravais_TYP *group,
                                        matrix_TYP *gram,
				       int *divisors,
				       ZZ_data_t *data,
				       ZZ_tree_t *tree,
				       int *projections,
				       int konst_flag);
extern matrix_TYP *ZZ_fget_data(ZZ_data_t * data,
					       ZZ_tree_t * tree,
					       char *file_name);
extern void ZZ_free_data(ZZ_data_t * data);
extern void ZZ_put_data(bravais_TYP *, ZZ_data_t * data, 
				       ZZ_tree_t * tree);
extern void ZZ_fput_data(ZZ_data_t * data, ZZ_tree_t * tree, 
					int ABBRUCH);
extern matrix_TYP *ZZ_p_vsolve(int *anz, matrix_TYP ** L_mat, 
					      matrix_TYP ** R_mat);
extern ZZ_node_t *ZZ_center(ZZ_data_t * data, 
					   ZZ_node_t * father,
					   int ii, int jj);
extern int ZZ_epimorphs(ZZ_data_t * data, int ii, int jj);
extern boolean ZZ_successor(ZZ_data_t * data, 
					   ZZ_node_t ** act);
extern int ZZ_ins_node(matrix_TYP * Gram,
				      ZZ_data_t * data, ZZ_tree_t * tree,
				      ZZ_node_t * father, ZZ_node_t * newnode,
				      int ii, int jj,
				      QtoZ_TYP *inzidenz,
				      int *nr,
				      int *NEU,
				      int *flagge,
				      int *g,
				      ZZ_node_t **nnn);
extern void ZZ_pick_epi(ZZ_data_t * data, int number, 
				       int ii, int jj);

#endif
