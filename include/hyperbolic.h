#ifdef __cplusplus
extern "C" {
#endif


#ifndef _HYPERBOLIC_H_
#define _HYPERBOLIC_H_

#ifndef _CARAT_TYPEDEF_H_
#include"typedef.h"
#endif

#ifdef __STDC__
/*-------------------------------------------------------------*\
| FILE: hyp_isom.c
\*-------------------------------------------------------------*/
extern matrix_TYP *hyperbolic_isometry(matrix_TYP *x1, matrix_TYP *x2,
     matrix_TYP *S);

/*-------------------------------------------------------------*\
| FILE: hyp_stabilizer.c
\*-------------------------------------------------------------*/
extern bravais_TYP *hyperbolic_stabilizer(matrix_TYP *x, matrix_TYP *S);
#else
/*-------------------------------------------------------------*\
| FILE: hyp_isom.c
\*-------------------------------------------------------------*/
extern matrix_TYP *hyperbolic_isometry();

/*-------------------------------------------------------------*\
| FILE: hyp_stabilizer.c
\*-------------------------------------------------------------*/
extern bravais_TYP *hyperbolic_stabilizer();

#endif
#endif


#ifdef __cplusplus
}
#endif

