#ifdef __cplusplus
extern "C" {
#endif



#ifndef _CONTRIB_H_
#define _CONTRIB_H_

#ifndef _CARAT_TYPEDEF_H_
#include"typedef.h"
#endif

#ifdef __STDC__

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

#else

/*-------------------------------------------------------------*\
| FILE: torsionfree.c
\*-------------------------------------------------------------*/

int *torsionfree();

/***************************************************************
| FILE: suche_kand.c
****************************************************************/

matrix_TYP *suche_kand ();

#endif
#endif


#ifdef __cplusplus
}
#endif

