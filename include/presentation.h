#ifndef _PRESENTATION_h_ 
#define _PRESENTATION_h_
#include "baum.h"
#include "tietzetrans.h"

#ifdef __STDC__
/**********************************************************************
|  FILE: neugram.c
***********************************************************************/
extern matrix_TYP *neugram(matrix_TYP *,
                           matrix_TYP *);

/**********************************************************************
|  FILE: KKGV.c
***********************************************************************/
extern int KKGV(int *,
                int);

/**********************************************************************
|  FILE: neukandidat.c
***********************************************************************/
extern matrix_TYP *neukandidat_vectors(matrix_TYP *,
                                       matrix_TYP *);

/**********************************************************************
|  FILE: modshort.c
***********************************************************************/
extern matrix_TYP *neumodshort_vectors( matrix_TYP *,
                                        matrix_TYP *,
                                        int );

/**********************************************************************
|  FILE: sscal.c
***********************************************************************/
extern void sscal_mul(matrix_TYP *mat,
                      int v);

/**********************************************************************
|  FILE: modshort.c
***********************************************************************/
extern matrix_TYP *modshort_vectors(matrix_TYP *,
                                    matrix_TYP *,
                                    int);

/**********************************************************************
|  FILE: T_polyeder.c
***********************************************************************/
extern polyeder_TYP *T_polyeder(matrix_TYP *,
                                matrix_TYP *,
                                matrix_TYP *);

/**********************************************************************
|  FILE: T_ungleich_wand.c
***********************************************************************/
extern wall_TYP **T_ungleich_wand(matrix_TYP *,
                                  matrix_TYP *,
                                  matrix_TYP *,
                                  int *);

/***************************************************************\
|  FILE : tietze_tools.c
\***************************************************************/
void put_presentation(presentation_TYP *pres,
                      char *file,
                      char *option);

void free_presentation(presentation_TYP *pres);


/***************************************************************\
|  FILE : sscal_mul.c
\***************************************************************/
void sscal_mul(matrix_TYP *mat,
               int v);

/***************************************************************\
|  FILE : fub_tools.c
\***************************************************************/
void Stab_word_mul(Stab_word_TYP *sw1,
                   Stab_word_TYP *sw2,
                   matrix_TYP **Mat,
                   matrix_TYP **Mat_inv,
                   int matanz,
                   Stab_word_TYP *erg);

void Stab_word_mul(Stab_word_TYP *sw1,
                   Stab_word_TYP *sw2,
                   matrix_TYP **Mat,
                   matrix_TYP **Mat_inv,
                   int matanz,
                   Stab_word_TYP *erg);

void free_Stab_word(Stab_word_TYP *SW);

/***************************************************************\
|  FILE: presentation_point_grp.c
\***************************************************************/
presentation_TYP *presentation_point_grp(bravais_TYP *G);

/***************************************************************\
|  FILE: barzen.c
\***************************************************************/
int *make_prod(word_TYP *word);

matrix_TYP *barzen(polyeder_TYP *Pol);

presentation_TYP* make_pres(anne_presentation_TYP *PR);

int* generate_Liste(polyeder_TYP *Pol,
                    int **Liste_inv);

/***************************************************************\
| FILE: fub_tools.c
\***************************************************************/
extern matrix_TYP *word_to_mat(int *word,
                               matrix_TYP **Mat,
                               matrix_TYP **Mat_inv,
                               int matanz);

extern  matrix_TYP *orb_word_to_mat(int *word,
                                    matrix_TYP **Mat,
                                    matrix_TYP **Mat_inv,
                                    int matanz,
                                    matrix_TYP *trans);

extern matrix_TYP *Stab_word_to_mat(Stab_word_TYP SW,
                                    matrix_TYP **Mat, 
                                    matrix_TYP **Mat_inv,
                                    int matanz,
                                    matrix_TYP *Tmat);

extern int *wand_ungl(matrix_TYP *M,
                      matrix_TYP *B,
                      matrix_TYP *Form);

extern int wand(int erz_anz,
                matrix_TYP **Erz,
                matrix_TYP **Erz_inv,
                ele_TYP ** orbit,
                int length,
                polyeder_TYP * P,
                matrix_TYP * Form,
                int **rel_kand);


extern void P_ausgabe(polyeder_TYP *P);

extern void free_ele(ele_TYP *orb); 

extern void Word_mul( word_TYP *sw1,
                      word_TYP *sw2,
                      matrix_TYP **Mat,
                      matrix_TYP **Mat_inv,
                      int matanz,
                      word_TYP *erg);


#else

/**********************************************************************
|  FILE: neugram.c
***********************************************************************/
extern matrix_TYP *neugram();

/**********************************************************************
|  FILE: KKGV.c
***********************************************************************/
extern int KKGV();

/**********************************************************************
|  FILE: neukandidat.c
***********************************************************************/
extern matrix_TYP *neukandidat_vectors();

/**********************************************************************
|  FILE: modshort.c
***********************************************************************/
extern matrix_TYP *neumodshort_vectors();

/**********************************************************************
|  FILE: sscal.c
***********************************************************************/
extern void sscal_mul();

/**********************************************************************
|  FILE: modshort.c
***********************************************************************/
extern matrix_TYP *modshort_vectors();

/**********************************************************************
|  FILE: T_polyeder.c
***********************************************************************/
extern polyeder_TYP *T_polyeder();

/**********************************************************************
|  FILE: T_ungleich_wand.c
***********************************************************************/
extern wall_TYP **T_ungleich_wand();

/***************************************************************\
|  FILE : tietze_tools.c
\***************************************************************/
void put_presentation();

void free_presentation();


/***************************************************************\
|  FILE : sscal_mul.c
\***************************************************************/
void sscal_mul();

/***************************************************************\
|  FILE : fub_tools.c
\***************************************************************/
void Stab_word_mul();

void Stab_word_mul();

void free_Stab_word();

/***************************************************************\
|  FILE: presentation_point_grp.c
\***************************************************************/
presentation_TYP *presentation_point_grp();

/***************************************************************\
|  FILE: barzen.c
\***************************************************************/
int *make_prod();

matrix_TYP *barzen();

presentation_TYP* make_pres();

int* generate_Liste();

/***************************************************************\
| FILE: fub_tools.c
\***************************************************************/
extern matrix_TYP *word_to_mat();

extern  matrix_TYP *orb_word_to_mat();

extern matrix_TYP *Stab_word_to_mat();

extern int *wand_ungl();

extern int wand();


extern void P_ausgabe();

extern void free_ele(); 

extern void Word_mul();

#endif
#endif

