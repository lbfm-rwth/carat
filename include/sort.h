#ifndef CARAT_SORT_H
#define CARAT_SORT_H

#include "typedef.h"

/*-------------------------------------------------------------*\
| FILE: compare.c 
\*-------------------------------------------------------------*/
extern int mat_comp(const matrix_TYP *m1, const matrix_TYP *m2);
extern int mat_col_comp(const matrix_TYP *m1, const matrix_TYP *m2);
extern int lower_triangular_mat_comp(const matrix_TYP *m1, const matrix_TYP *m2);
extern int vec_comp(const int *v1, const int *v2, int dim);
extern int pointer_mat_comp(const int * const *m1, const int * const *m2, int rows, int cols);
extern int pointer_lower_triangular_mat_comp(const int * const *m1, const int * const *m2, int n, int m);

/*-------------------------------------------------------------*\
| FILE: quicksort.c 
\*-------------------------------------------------------------*/
extern void mat_quicksort(matrix_TYP **M, int inf, int sup, int (*comp)(const matrix_TYP *, const matrix_TYP *));
extern void vec_quicksort(int **v, int inf, int sup, int dim, int (*comp)(const int *, const int *, int));
extern void pointer_mat_quicksort(int ***v, int inf, int sup, int rows,
     int cols, int (*comp)(const int * const *, const int * const *, int, int));

/*-------------------------------------------------------------*\
| FILE: search.c 
\*-------------------------------------------------------------*/
extern int mat_search(const matrix_TYP *M, const matrix_TYP * const *List, int List_no,
     int (*comp)(const matrix_TYP *, const matrix_TYP *));
extern int vec_search(const matrix_TYP *M, const matrix_TYP * const *List, int List_no,
     int dim, int (*comp)(const int *, const int *, int));
extern int pointer_mat_search(const matrix_TYP * const *M, const matrix_TYP * const * const  *List, int List_no,
     int rows, int cols, int (*comp)(const int **, const int **, int, int));

#endif
