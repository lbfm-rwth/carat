#ifdef __cplusplus
extern "C" {
#endif

#ifndef _PRESENTATION_h_ 
#define _PRESENTATION_h_

/**********************************************************************
|  FILE: presentation.c
***********************************************************************/
extern matrix_TYP *pres(bahn **s,
                        bravais_TYP *G,
                        int *OPT);

extern void normalize_word(int *w);

/**********************************************************************
|  FILE: mapped_word.c
***********************************************************************/
matrix_TYP *mapped_word(int *w,
                        matrix_TYP **A,
                        matrix_TYP **AINV);

/**********************************************************************
|  FILE: check_base.c
***********************************************************************/
void check_base(bahn **s,
                bravais_TYP *G);

/**********************************************************************
|  FILE: put_word.c
***********************************************************************/
void put_word(int *w,
              char *O);

#endif

#ifdef __cplusplus
}
#endif
