#ifdef __cplusplus
extern "C" {
#endif


#ifndef _SYMM_H_
#define _SYMM_H_

#ifndef _CARAT_TYPEDEF_H_
#include"typedef.h"
#endif

#ifdef __STDC__
/*-------------------------------------------------------------*\
| FILE: dsylv.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *dsylv(matrix_TYP *M);
extern int definite_test(matrix_TYP *M);

/*-------------------------------------------------------------*\
| FILE: rest_short.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *rest_short(matrix_TYP *mat, int *restvec, int rkgv,
     int zaehler, int nenner, int find_opt, int count_opt, int *anz);

/*-------------------------------------------------------------*\
| FILE: short.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *short_vectors(matrix_TYP *mat, int length,
     int lengthmin, int find_opt, int count_opt, int *anz);

/*-------------------------------------------------------------*\
| FILE: shortest.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *shortest(matrix_TYP *mat, int *min_norm);

#else
/*-------------------------------------------------------------*\
| FILE: dsylv.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *dsylv();
extern int definite_test();

/*-------------------------------------------------------------*\
| FILE: rest_short.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *rest_short();

/*-------------------------------------------------------------*\
| FILE: short.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *short_vectors();

/*-------------------------------------------------------------*\
| FILE: shortest.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *shortest();

#endif
#endif

#ifdef __cplusplus
}
#endif

