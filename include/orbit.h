#ifndef _ORBIT_H_
#define _ORBIT_H_

#ifndef _CARAT_TYPEDEF_H_
#include"typedef.h"
#endif

#ifdef __STDC__
/*-------------------------------------------------------------*\
| FILE: orb_division.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *orbit_representatives(matrix_TYP **M, int Manz,
     bravais_TYP *G, int *option, int *orbit_no, int is_sorted);

/*-------------------------------------------------------------*\
| FILE: orb_alg.c 
\*-------------------------------------------------------------*/
extern int *make_orbit_options();
extern matrix_TYP **orbit_alg(matrix_TYP *M, bravais_TYP *G, bravais_TYP *S,
     int *option, int *length);
extern void free_baum(struct baum *p);

/*-------------------------------------------------------------*\
| FILE: orbit_subdivision.c 
\*-------------------------------------------------------------*/
extern int *orbit_subdivision(matrix_TYP *vecs, bravais_TYP *G, int *orbit_no);

/*-------------------------------------------------------------*\
| FILE: vec_orbit_division.c 
\*-------------------------------------------------------------*/
extern int *vec_orbit_division( int **V, int V_no, bravais_TYP *G,
     int orbit_no);

/*-------------------------------------------------------------*\
| FILE: row_spin.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *row_spin(matrix_TYP *x,matrix_TYP **G,int no,int option);
extern bravais_TYP *representation_on_lattice(matrix_TYP *x,bravais_TYP *G,
                                       int option);
matrix_TYP *translation_lattice(matrix_TYP **G,int number,matrix_TYP *P);

#else
/*-------------------------------------------------------------*\
| FILE: orb_division.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *orbit_representatives();

/*-------------------------------------------------------------*\
| FILE: orb_alg.c 
\*-------------------------------------------------------------*/
extern int *make_orbit_options();
extern matrix_TYP **orbit_alg();
extern void free_baum();

/*-------------------------------------------------------------*\
| FILE: orbit_subdivision.c 
\*-------------------------------------------------------------*/
extern int *orbit_subdivision();

/*-------------------------------------------------------------*\
| FILE: vec_orbit_division.c 
\*-------------------------------------------------------------*/
extern int *vec_orbit_division();

/*-------------------------------------------------------------*\
| FILE: row_spin.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *row_spin();
extern bravais_TYP *representation_on_lattice();
matrix_TYP *translation_lattice();

#endif
#endif
