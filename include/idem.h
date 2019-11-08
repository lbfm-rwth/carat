#ifdef __cplusplus
extern "C" {
#endif

#ifndef _IDEM_H_
#define _IDEM_H_

#include "typedef.h"

/*************************************************************************
*  FILE : almost_decomposable_lattice.c
**************************************************************************/

matrix_TYP *almost_decomposable_lattice(bravais_TYP *G);

/*************************************************************************
*  FILE : centr.c
**************************************************************************/

matrix_TYP **idempotente(matrix_TYP **gen,int gen_no,matrix_TYP *form,
             int *anz,int *dimc,int *dimcc);

matrix_TYP **solve_endo(matrix_TYP **A,matrix_TYP **B,int anz,int *dim);


matrix_TYP *zeros(matrix_TYP *A);

/*************************************************************************
* FILE: min_pol.c
**************************************************************************/

matrix_TYP *min_pol(matrix_TYP *A);

/*************************************************************************
*  FILE: symbol.c
**************************************************************************/

typedef struct { bravais_TYP *group;
                 matrix_TYP **centralizer;
                 int dimc;
                 matrix_TYP **ccentralizer;
                 int dimcc;
                 matrix_TYP *lattice;
               } constituent;

char *symbol(bravais_TYP *G,matrix_TYP *F);

/*************************************************************************
*  FILE: bravais_catalog.c
**************************************************************************/

symbol_out *read_symbol_from_string(const char *symb);

bravais_TYP *catalog_number(bravais_TYP *G,const char *symb,matrix_TYP **TR,
                            int *almost,int *zclass);

/*************************************************************************
*  FILE: v4_catalog.c
**************************************************************************/

bravais_TYP *catalog_number_v4(bravais_TYP *G,const char *symb,matrix_TYP **TR,
                               int *almost,int *zclass);

/*************************************************************************
*  FILE: z_equivalent.c
**************************************************************************/

matrix_TYP *z_equivalent(bravais_TYP *G,
                         bravais_TYP **G_tr,
                         bravais_TYP *H);

#endif

#ifdef __cplusplus
}
#endif
