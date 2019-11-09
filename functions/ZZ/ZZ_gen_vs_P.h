#ifndef CARAT_ZZ_GEN_VS_P_H
#define CARAT_ZZ_GEN_VS_P_H

#include "typedef.h"
#include "matrix.h"
#include "ZZ_P.h"

extern matrix_TYP **ZZ_mod_gen(matrix_TYP ** gen, int num);
extern matrix_TYP *ZZ_vec_bahn(int *vec, matrix_TYP ** gen, 
					      int num);
extern matrix_TYP *ZZ_gen_vs(int prim, int dim);

#endif
