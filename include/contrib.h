#ifdef __cplusplus
extern "C" {
#endif

#ifndef _CONTRIB_H_
#define _CONTRIB_H_

#include "typedef.h"

/*-------------------------------------------------------------*\
| FILE: torsionfree.c
\*-------------------------------------------------------------*/

int *torsionfree(bravais_TYP *R,
                 int *order_out,
                 int *number_of_conjugacy_classes);

/***************************************************************
| FILE: suche_kand.c
****************************************************************/

matrix_TYP *suche_kand (bravais_TYP *Gen_A, bravais_TYP *Gen_B);

#endif

#ifdef __cplusplus
}
#endif
