#if !defined(_ZZ_H) && !defined(_ZZ_IRR_CONST_H)
#define _ZZ_IRR_CONST_H

#include "matrix.h"
#include "ZZ_P.h"

extern matrix_TYP **ZZ_irr_const(matrix_TYP ** generators, 
						int num_gens,
						int p, int *num_irr_const);

#endif
