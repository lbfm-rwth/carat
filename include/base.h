#ifdef __cplusplus
extern "C" {
#endif



#include"gmp.h"

#ifndef _BASE_H_
#define _BASE_H_

#ifdef __STDC__

/*************************************************************************
@ FILE: base.c
**************************************************************************/
extern void free_tree(struct tree *p);

extern int hash_mat(struct tree *p,matrix_TYP **v,matrix_TYP *x,int pos);

extern void init_bahn(bahn *a);

extern void free_bahn(bahn *a);

extern bahn **strong_generators(matrix_TYP **base,bravais_TYP *U,int OPT);

extern matrix_TYP **get_base(bravais_TYP *U);

extern int is_element(matrix_TYP *x,bravais_TYP *G,bahn **strong,int **w);

extern matrix_TYP **normalizer_in_N(bravais_TYP *U,bravais_TYP *N,int *anz,
                                    int finite_flag);

extern int size(bahn **a);

extern void extend_bahn(bahn **a);

extern matrix_TYP **get_base(bravais_TYP *U);

extern int red_gen(bravais_TYP *G,matrix_TYP **base,bahn ***strong,int i);

matrix_TYP *conjugated(bravais_TYP *G,bravais_TYP *H,
                       matrix_TYP **N,int Nanz,bahn **strong);

/*************************************************************************
@ FILE: base2.c
**************************************************************************/
int strong_generators_2(matrix_TYP **base,bravais_TYP *U,matrix_TYP ***K,
                        int *anz,MP_INT *mp);

#else

/*************************************************************************
@ FILE: base.c
**************************************************************************/
extern void free_tree();

extern int hash_mat();

extern void init_bahn();

extern void free_bahn();

extern bahn **strong_generators();

extern matrix_TYP **get_base();

extern int is_element();

extern matrix_TYP **normalizer_in_N();

extern int size();

extern void extend_bahn();

extern matrix_TYP **get_base();

extern int red_gen();

matrix_TYP *conjugated();

/*************************************************************************
@ FILE: base2.c
**************************************************************************/
int strong_generators_2();

#endif
#endif

#ifdef __cplusplus
}
#endif


