#ifndef _ZASSEN_H_
#define _ZASSEN_H_

#ifndef _CARAT_TYPEDEF_H_
#include"typedef.h"
#endif

int INFO_LEVEL;

typedef struct {
       long last;
       long speicher;
       long faktor;
       long *pointer;} word;


#ifdef __STDC__
extern void matrix_2_word(matrix_TYP *matrix,word *relator,long zeile);

extern matrix_TYP **cohomology(long *dim,matrix_TYP **mat,matrix_TYP **matinv,
                word *relator,int erzeuger,int relatoren);

extern int wordfree(word *a);

matrix_TYP *scalar(long n, long a);


/**********************************************************************\
| FILE: cobundary.c
\**********************************************************************/

void coboundary(bravais_TYP *G,
                matrix_TYP *C,
                matrix_TYP *T);

/**********************************************************************\
| FILE: cong_solve.c
\**********************************************************************/

matrix_TYP **cong_solve(matrix_TYP *A);

/**********************************************************************\
| FILE: convert_cocycle_to_column.c
\**********************************************************************/

void convert_cocycle_to_column(matrix_TYP **Y,
                               int number,
                               int dim,
                               int gen_no);

/**********************************************************************\
| FILE: convert_to_cozycle.c
\**********************************************************************/

extern matrix_TYP *convert_to_cozycle(matrix_TYP *x,
                                      matrix_TYP *cozycle,
                                      matrix_TYP *D);

/**********************************************************************\
| FILE: put_cocycle.c
\**********************************************************************/

void put_cocycle(matrix_TYP *COZ,
                 int dim,
                 int number,
                 char *file,
                 char *comment);

/*********************************************************************\
| FILE: normalop.c
\*********************************************************************/

extern matrix_TYP *normalop(matrix_TYP *cozycle,
                            matrix_TYP *D,
                            matrix_TYP *R,
                            bravais_TYP *G,
                            matrix_TYP *N,
                            int opt);

extern matrix_TYP **extensions(matrix_TYP *cozycle,
                               matrix_TYP *D,
                               matrix_TYP *R,
                               bravais_TYP *G,
                               int **lengths,
                               MP_INT **names,
                               int *number_of_orbits,
                               int option);

extern void no_of_extensions(matrix_TYP *cozycle,matrix_TYP *D,
                        matrix_TYP *R,bravais_TYP *G,MP_INT *no);

extern matrix_TYP **identify(matrix_TYP *cozycle,
                             matrix_TYP *D,
                             matrix_TYP *R,
                             bravais_TYP *G,
                             matrix_TYP **extension,
                             MP_INT *a,int number,
                             int transform_flag,
                             int ***WORDS,
                             int *NUMBER_OF_WORDS);


/**********************************************************************\
| FILE: reverse_valuation.c
\**********************************************************************/

extern matrix_TYP *reverse_valuation(MP_INT *val,matrix_TYP *D);

/**********************************************************************\
| FILE: reget_gen.c
\**********************************************************************/

extern matrix_TYP *reget_gen(matrix_TYP **map,int number,bravais_TYP *G,
                             int **words, int word_flag);
#else

extern void matrix_2_word();
extern matrix_TYP **cohomology();
extern int wordfree();
matrix_TYP *scalar();

/**********************************************************************\
| FILE: cobundary.c
\**********************************************************************/

void coboundary();

/**********************************************************************\
| FILE: cong_solve.c
\**********************************************************************/

matrix_TYP **cong_solve();


/**********************************************************************\
| FILE: convert_cocycle_to_column.c
\**********************************************************************/

void convert_cocycle_to_column();

/**********************************************************************\
| FILE: convert_to_cozycle.c
\**********************************************************************/

extern matrix_TYP *convert_to_cozycle();


/**********************************************************************\
| FILE: put_cocycle.c
\**********************************************************************/

void put_cocycle();

/*********************************************************************\
| FILE: normalop.c
\*********************************************************************/

extern matrix_TYP *normalop();

extern matrix_TYP **extensions();

extern void no_of_extensions();

extern matrix_TYP **identify();

/**********************************************************************\
| FILE: reverse_valuation.c
\**********************************************************************/

extern matrix_TYP *reverse_valuation();

/**********************************************************************\
| FILE: reget_gen.c
\**********************************************************************/

extern matrix_TYP *reget_gen();

#endif
#endif
