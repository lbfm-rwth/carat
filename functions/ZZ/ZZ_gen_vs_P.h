#if !defined(_ZZ_H) && !defined(_ZZ_GEN_VS_P_H)
#define ZZ_GEN_VS_P_H

#include "typedef.h"
#include "matrix.h"
#include "ZZ_P.h"

extern matrix_TYP **ZZ_mod_gen _ZZ_P_PROTO_( (matrix_TYP ** gen, int num) );
extern matrix_TYP *ZZ_vec_bahn _ZZ_P_PROTO_( (int *vec, matrix_TYP ** gen, 
					      int num) );
extern matrix_TYP *ZZ_gen_vs _ZZ_P_PROTO_( (int prim, int dim) );

#endif
