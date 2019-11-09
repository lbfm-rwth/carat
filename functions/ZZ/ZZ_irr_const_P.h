#ifndef CARAT_ZZ_IRR_CONST_H
#define CARAT_ZZ_IRR_CONST_H

#include "matrix.h"
#include "ZZ_P.h"

extern matrix_TYP **ZZ_irr_const(matrix_TYP ** generators, 
						int num_gens,
						int p, int *num_irr_const);

#endif
