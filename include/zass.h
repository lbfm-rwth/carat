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

/*********************************************************************\
| FILE: normalop.c
\*********************************************************************/

extern matrix_TYP *normalop(matrix_TYP *cozycle,matrix_TYP *D,matrix_TYP *R,
                     bravais_TYP *G,matrix_TYP *N);

extern matrix_TYP **extensions(matrix_TYP *cozycle,matrix_TYP *D,matrix_TYP *R,
             bravais_TYP *G,int **lengths,MP_INT **names,int *number_of_orbits);

extern void no_of_extensions(matrix_TYP *cozycle,matrix_TYP *D,
                        matrix_TYP *R,bravais_TYP *G,MP_INT *no);

extern matrix_TYP **identify(matrix_TYP *cozycle,matrix_TYP *D,matrix_TYP *R,
            bravais_TYP *G, matrix_TYP **extension,MP_INT *a,int number,
            int transform_flag);

/**********************************************************************\
| FILE: reget_gen.c
\**********************************************************************/

extern matrix_TYP *reget_gen(matrix_TYP **map,int number,bravais_TYP *G,
                             int **words, int word_flag);
#else

extern void matrix_2_word();
extern matrix_TYP **cohomology();
extern int wordfree();

/*********************************************************************\
| FILE: normalop.c
\*********************************************************************/

extern matrix_TYP *normalop();

extern matrix_TYP **extensions();

extern void no_of_extensions();

extern matrix_TYP **identify();

/**********************************************************************\
| FILE: reget_gen.c
\**********************************************************************/

extern matrix_TYP *reget_gen();

#endif
#endif
