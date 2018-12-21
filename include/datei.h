#ifdef __cplusplus
extern "C" {
#endif


#ifndef _DATEI_H_
#define _DATEI_H_

#ifndef _CARAT_TYPEDEF_H_
#include"typedef.h"
#endif

typedef struct{
    bravais_TYP *grp;
    char *symbol;
    int almost;
    int zclass;
    int alpha;
    int N_orbits;
    matrix_TYP **TR;
} lattice_element;

#ifdef __STDC__

/*-------------------------------------------------------------*\
| FILE: brav_from_datei.c
\*-------------------------------------------------------------*/
bravais_TYP *brav_from_datei(char *symb,int almost,int zclass);

/*-------------------------------------------------------------*\
| FILE: free_bravais.c 
\*-------------------------------------------------------------*/
extern void free_bravais( bravais_TYP *grp);

/*-------------------------------------------------------------*\
| FILE: get_symbol.c 
\*-------------------------------------------------------------*/
extern symbol_out *get_symbol ( char *file_name);

/*-------------------------------------------------------------*\
| FILE: get_zentr.c 
\*-------------------------------------------------------------*/
extern void get_zentr( symbol_out *B);

/*-------------------------------------------------------------*\
| FILE: gittstab.c 
\*-------------------------------------------------------------*/
extern bravais_TYP *gittstab( bravais_TYP *grp, matrix_TYP *X);

extern bravais_TYP *Z_class( bravais_TYP *B, matrix_TYP *zen);

/*-------------------------------------------------------------*\
| FILE: gittstabneu.c
\*-------------------------------------------------------------*/
extern bravais_TYP *gittstabneu( bravais_TYP *grp, matrix_TYP *X);

/*-------------------------------------------------------------*\
| FILE: read_symbol.c 
\*-------------------------------------------------------------*/
extern symbol_out *read_symbol(char *file_name);

/*-------------------------------------------------------------*\
| FILE: right_order.c
\*-------------------------------------------------------------*/
extern void right_order(char *string);

/*-------------------------------------------------------------*\
| FILE: lattice_tools.c
\*-------------------------------------------------------------*/

lattice_element *init_lattice_element();
void free_lattice_element(lattice_element *x);
lattice_element *fget_lattice_element(FILE *F,int OPTION);
void fput_lattice_element(lattice_element *E,FILE *F);

/*-------------------------------------------------------------*\
| FILE: lattice.c
\*-------------------------------------------------------------*/

lattice_element **lattice(char *symb,int dim,int almost,int zclass,int *no,
                          int OPTION);

/*------------------------------------------------------------*\
| FILE: super_lattice.c
\*------------------------------------------------------------*/

lattice_element **super_lattice(char *symb,int dim,int almost,int zclass,
		int *no, int OPTION);

/*------------------------------------------------------------*\
| FILE: get_data_dir.c
\*------------------------------------------------------------*/

void get_data_dir(char *result, const char *str);

#else

/*-------------------------------------------------------------*\
| FILE: brav_from_datei.c
\*-------------------------------------------------------------*/
bravais_TYP *brav_from_datei();

/*-------------------------------------------------------------*\
| FILE: free_bravais.c 
\*-------------------------------------------------------------*/
extern void free_bravais();

/*-------------------------------------------------------------*\
| FILE: get_symbol.c 
\*-------------------------------------------------------------*/
extern symbol_out *get_symbol();

/*-------------------------------------------------------------*\
| FILE: get_zentr.c 
\*-------------------------------------------------------------*/
extern void get_zentr();

/*-------------------------------------------------------------*\
| FILE: gittstab.c 
\*-------------------------------------------------------------*/
extern bravais_TYP *gittstab();

extern bravais_TYP *Z_class();

/*-------------------------------------------------------------*\
| FILE: gittstabneu.c
\*-------------------------------------------------------------*/
extern bravais_TYP *gittstabneu();

/*-------------------------------------------------------------*\
| FILE: read_symbol.c 
\*-------------------------------------------------------------*/
symbol_out *read_symbol();

/*-------------------------------------------------------------*\
| FILE: right_order.c
\*-------------------------------------------------------------*/
extern void right_order();

/*--------------------------------------------------------------------------*\
| FILE: lattice_tools.c
\*--------------------------------------------------------------------------*/

lattice_element *init_lattice_element();
void free_lattice_element();
lattice_element *fget_lattice_element();
void fput_lattice_element();

/*-------------------------------------------------------------*\
| FILE: lattice.c
\*-------------------------------------------------------------*/

lattice_element **lattice(char *symb,int dim,int almost,int zclass,int *no,
                          int OPTION);

/*------------------------------------------------------------*\
| FILE: get_data_dir.c
\*------------------------------------------------------------*/

void get_data_dir();

#endif
#endif


#ifdef __cplusplus
}
#endif

