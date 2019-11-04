#ifdef __cplusplus
extern "C" {
#endif


#ifndef _SORT_H_
#define _SORT_H_

#include "typedef.h"

/*-------------------------------------------------------------*\
| FILE: compare.c 
\*-------------------------------------------------------------*/
extern int mat_comp(matrix_TYP *m1, matrix_TYP *m2);
extern int mat_col_comp(matrix_TYP *m1, matrix_TYP *m2);
extern int lower_triangular_mat_comp(matrix_TYP *m1,matrix_TYP *m2);
extern int vec_comp(int *v1, int *v2, int dim);
extern int pointer_mat_comp(int **m1, int **m2, int rows, int cols);
extern int pointer_lower_triangular_mat_comp(int **m1,int **m2, int n, int m);

/*-------------------------------------------------------------*\
| FILE: quicksort.c 
\*-------------------------------------------------------------*/
extern void mat_quicksort(matrix_TYP **M, int inf, int sup, int (*comp)());
extern void vec_quicksort(int **v, int inf, int sup, int dim, int (*comp)());
extern void pointer_mat_quicksort(int ***v, int inf, int sup, int rows,
     int cols, int (*comp)());

/*-------------------------------------------------------------*\
| FILE: search.c 
\*-------------------------------------------------------------*/
extern int mat_search(matrix_TYP *M, matrix_TYP **List, int List_no,
     int (*comp)());
extern int vec_search(matrix_TYP *M, matrix_TYP **List, int List_no,
     int dim, int (*comp)());
extern int pointer_mat_search(matrix_TYP **M, matrix_TYP ***List, int List_no,
     int rows, int cols, int (*comp)());

#endif

#ifdef __cplusplus
}
#endif
