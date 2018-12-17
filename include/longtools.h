#ifdef __cplusplus
extern "C" {
#endif


#ifndef _LONGTOOLS_H_
#define _LONGTOOLS_H_

#ifndef _CARAT_TYPEDEF_H_
#include"typedef.h"
#endif

#ifndef _GMP_H_
#include "gmp.h"
#endif

#ifdef __STDC__
/*-------------------------------------------------------------*\
| FILE: MP_conv_mat.c
\*-------------------------------------------------------------*/
extern MP_INT **matrix_to_MP_mat(matrix_TYP *M);
extern matrix_TYP *MP_mat_to_matrix(MP_INT **M, int rows, int cols);
extern void write_MP_mat_to_matrix(matrix_TYP *Mat, MP_INT **mp);
extern MP_INT **init_MP_mat(int rows, int cols);
extern void free_MP_mat(MP_INT **M, int rows, int cols);

/*-------------------------------------------------------------*\
| FILE: MP_gauss.c
\*-------------------------------------------------------------*/
extern int MP_trf_gauss(MP_INT **M, MP_INT **Trf, int rows, int cols);
extern int MP_row_gauss(MP_INT **M, int rows, int cols);
extern int MP_row_gauss_simultaneous(MP_INT **M, int rows, int cols,
     MP_INT **B, int Bcols);
extern void MP_row_gauss_reverse(MP_INT **A,int rows,int cols,int option);

/*-------------------------------------------------------------*\
| FILE: MP_hnf.c
\*-------------------------------------------------------------*/
extern int MP_trf_hnf(MP_INT **M, MP_INT **Trf, int rows, int cols);
extern int MP_hnf(MP_INT **M, int rows, int cols);
extern int MP_hnf_simultaneous(MP_INT **M, int rows, int cols, MP_INT **B,
     int Bcols);

/*-------------------------------------------------------------*\
| FILE: MP_pair_red.c
\*-------------------------------------------------------------*/
extern void MP_pair_red(MP_INT **G, MP_INT **T, int n);

/*-------------------------------------------------------------*\
| FILE: MP_red_sort.c
\*-------------------------------------------------------------*/
extern void MP_reduction_sort(MP_INT **G,MP_INT **T,int n);

/*-------------------------------------------------------------*\
| FILE: MP_solve.c
\*-------------------------------------------------------------*/
extern MP_INT ***MP_solve_mat(MP_INT **M, int rows, int cols, MP_INT **B,
     int Bcols, int *X1cols, MP_INT *X0kgv);

/*-------------------------------------------------------------*\
| FILE: long_elt.c
\*-------------------------------------------------------------*/
extern matrix_TYP *long_elt_mat(matrix_TYP *left_trans,
                                matrix_TYP *Mat,
                                matrix_TYP *right_trans);

/*-------------------------------------------------------------------*\
|  FILE: long_gauss.c
\*-------------------------------------------------------------------*/
extern int long_row_gauss(matrix_TYP *Mat);
extern int long_row_basis(matrix_TYP *Mat,int flag);
extern int long_row_trf_gauss(matrix_TYP *M, matrix_TYP *T);
extern int long_row_gauss_simultaneous(matrix_TYP *A, matrix_TYP *B);

/*-------------------------------------------------------------------*\
|  FILE: long_hnf.c
\*-------------------------------------------------------------------*/
extern int long_row_hnf(matrix_TYP *Mat);
extern int long_col_hnf(matrix_TYP *Mat);
extern int long_row_trf_hnf(matrix_TYP *M, matrix_TYP *T);
extern int long_row_hnf_simultaneous(matrix_TYP *A, matrix_TYP *B);

/*-------------------------------------------------------------*\
| FILE: long_kernel_mat.c
\*-------------------------------------------------------------*/
extern matrix_TYP *long_kernel_mat(matrix_TYP *A);

/*-------------------------------------------------------------*\
| FILE: long_mat_inv.c
\*-------------------------------------------------------------*/
extern matrix_TYP *long_mat_inv(matrix_TYP *A);

/*-------------------------------------------------------------*\
| FILE: long_qbase.c
\*-------------------------------------------------------------*/
extern matrix_TYP *long_qbase();

/*-------------------------------------------------------------*\
| FILE: long_rein_mat.c
\*-------------------------------------------------------------*/
extern matrix_TYP *long_rein_mat(matrix_TYP *M);
extern int long_rein_formspace(matrix_TYP **forms,int number,int option);

/*-------------------------------------------------------------*\
| FILE: long_solve_mat.c
\*-------------------------------------------------------------*/
extern matrix_TYP **long_solve_mat(matrix_TYP *A, matrix_TYP *B);

/*-------------------------------------------------------------*\
| FILE: dump_MP_mat.c
\*-------------------------------------------------------------*/
extern void dump_MP_mat(MP_INT **Mat, int rows, int cols, char *comment);
  
#else
/*-------------------------------------------------------------*\
| FILE: MP_conv_mat.c
\*-------------------------------------------------------------*/
extern MP_INT **matrix_to_MP_mat();
extern matrix_TYP *MP_mat_to_matrix();
extern void write_MP_mat_to_matrix();
extern MP_INT **init_MP_mat();
extern void free_MP_mat();

/*-------------------------------------------------------------*\
| FILE: MP_gauss.c
\*-------------------------------------------------------------*/
extern int MP_trf_gauss();
extern int MP_row_gauss();
extern int MP_row_gauss_simultaneous();
extern void MP_row_gauss_reverse();

/*-------------------------------------------------------------*\
| FILE: MP_hnf.c
\*-------------------------------------------------------------*/
extern int MP_trf_hnf();
extern int MP_hnf();
extern int MP_hnf_simultaneous();

/*-------------------------------------------------------------*\
| FILE: MP_pair_red.c
\*-------------------------------------------------------------*/
extern void MP_pair_red();

/*-------------------------------------------------------------*\
| FILE: MP_red_sort.c
\*-------------------------------------------------------------*/
extern void MP_reduction_sort();

/*-------------------------------------------------------------*\
| FILE: MP_solve.c
\*-------------------------------------------------------------*/
extern MP_INT ***MP_solve_mat();

/*-------------------------------------------------------------*\
| FILE: long_elt.c
\*-------------------------------------------------------------*/
extern matrix_TYP *long_elt_mat();

/*-------------------------------------------------------------------*\
|  FILE: long_gauss.c
\*-------------------------------------------------------------------*/
extern int long_row_gauss();
extern int long_row_basis();
extern int long_row_trf_gauss();
extern int long_row_gauss_simultaneous();

/*-------------------------------------------------------------------*\
|  FILE: long_hnf.c
\*-------------------------------------------------------------------*/
extern int long_row_hnf();
extern int long_row_trf_hnf();
extern int long_row_hnf_simultaneous();

/*-------------------------------------------------------------*\
| FILE: long_kernel_mat.c
\*-------------------------------------------------------------*/
extern matrix_TYP *long_kernel_mat();

/*-------------------------------------------------------------*\
| FILE: long_mat_inv.c
\*-------------------------------------------------------------*/
extern matrix_TYP *long_mat_inv();

/*-------------------------------------------------------------*\
| FILE: long_qbase.c
\*-------------------------------------------------------------*/
extern matrix_TYP *long_qbase(matrix_TYP *Mat);

/*-------------------------------------------------------------*\
| FILE: long_rein_mat.c
\*-------------------------------------------------------------*/
extern matrix_TYP *long_rein_mat();
extern int long_rein_formspace();

/*-------------------------------------------------------------*\
| FILE: long_solve_mat.c
\*-------------------------------------------------------------*/
extern matrix_TYP **long_solve_mat();

/*-------------------------------------------------------------*\
| FILE: dump_MP_mat.c
\*-------------------------------------------------------------*/
extern void dump_MP_mat();

#endif
#endif


#ifdef __cplusplus
}
#endif

