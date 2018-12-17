#ifdef __cplusplus
extern "C" {
#endif


#ifndef _POLYEDER_H_
#define _POLYEDER_H_

#ifndef _CARAT_TYPEDEF_H_
#include"typedef.h"
#endif

#ifdef __STDC__
/*-------------------------------------------------------------*\
| FILE: first_polyeder.c 
\*-------------------------------------------------------------*/
extern polyeder_TYP *first_polyeder(wall_TYP **mauern, int anz);

/*-------------------------------------------------------------*\
| FILE: polyeder_tools.c 
\*-------------------------------------------------------------*/
extern vertex_TYP *init_vertex(int dim, int wall_no) ;
extern word_TYP *init_word(int dim);
extern wall_TYP *init_wall(int dim);
extern polyeder_TYP *init_polyeder(int vert_no, int wall_no);
extern polyeder_TYP *get_polyeder(char *file_name);
extern void put_polyeder(polyeder_TYP *F);
extern int wall_times_vertex(wall_TYP *w, vertex_TYP *v);
extern void free_vertex(vertex_TYP **v);
extern void free_word(word_TYP *word);
extern void free_wall(wall_TYP **v);
extern wall_TYP *mat_to_wall(matrix_TYP *M);
extern void normal_wall(wall_TYP *v);
extern void normal_vertex(vertex_TYP *v);
extern int is_vertex_of_wallno(vertex_TYP *v, int w);
extern word_TYP *copy_word(word_TYP *w);
extern wall_TYP *copy_wall(wall_TYP *w);
extern void free_polyeder(polyeder_TYP *P);

/*-------------------------------------------------------------*\
| FILE: polyeder_to_vecs.c 
\*-------------------------------------------------------------*/
extern matrix_TYP **polyeder_to_vecs(polyeder_TYP *P);

/*-------------------------------------------------------------*\
| FILE: refine_polyeder.c 
\*-------------------------------------------------------------*/
extern int refine_polyeder(polyeder_TYP *F, wall_TYP *h);

/*-------------------------------------------------------------*\
| FILE: fuber_tools.c 
\*-------------------------------------------------------------*/
void free_vertex_fuber(vertex_TYP**);
void free_wall_fuber(wall_TYP**);
int wall_times_vertex_fuber(wall_TYP*, vertex_TYP*);

#else
/*-------------------------------------------------------------*\
| FILE: first_polyeder.c 
\*-------------------------------------------------------------*/
extern polyeder_TYP *first_polyeder();

/*-------------------------------------------------------------*\
| FILE: polyeder_tools.c 
\*-------------------------------------------------------------*/
extern vertex_TYP *init_vertex();
extern wall_TYP *init_wall();
extern polyeder_TYP *init_polyeder();
extern polyeder_TYP *get_polyeder();
extern void put_polyeder();
extern int wall_times_vertex();
extern void free_vertex();
extern void free_word();
extern void free_wall();
extern wall_TYP *mat_to_wall();
extern void normal_wall();
extern void normal_vertex();
extern int is_vertex_of_wallno();
extern wall_TYP *copy_wall();
extern void free_polyeder();

/*-------------------------------------------------------------*\
| FILE: polyeder_to_vecs.c 
\*-------------------------------------------------------------*/
extern matrix_TYP **polyeder_to_vecs();

/*-------------------------------------------------------------*\
| FILE: refine_polyeder.c 
\*-------------------------------------------------------------*/
extern int refine_polyeder();

/*-------------------------------------------------------------*\
| FILE: fuber_tools.c 
\*-------------------------------------------------------------*/
void free_vertex_fuber();
void free_wall_fuber();
int wall_times_vertex_fuber();

#endif
#endif

#ifdef __cplusplus
}
#endif

