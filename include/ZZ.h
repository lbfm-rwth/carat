#ifndef _ZZ_H
#define _ZZ_H

#include <stdlib.h>		/* for atoi */
#include "typedef.h"
#include "tools.h"
#include "matrix.h"
#include "getput.h"
#include "datei.h"

#undef ZZ_PRIVATE

#ifdef __STDC__
extern void ZZ (bravais_TYP * group, matrix_TYP * gram, int *divisors, 
	       char *options, FILE* outputfile);

extern bravais_TYP **q2z(bravais_TYP *G,int *number,int ADFLAG);
#else
extern void ZZ ();
extern bravais_TYP **q2z();
#endif

extern int NUMBER;		/* Abbruch nach NUMBER Zentrierungen */
extern int LEVEL;		/* Abbruch nach Iterationszahl LEVEL */

#endif /* _ZZ_H */