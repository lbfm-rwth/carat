
/* check that this routine is included only once */
#ifndef _MALLOC_

#define _MALLOC_

#define m_alloc m_alloc_d1
#define c_alloc c_alloc_d1
#define re_alloc re_alloc_d1
#define fr_ee fr_ee_d1

/* decide which diagnostics we want */

#ifdef DIAG1
/* take this definition to have moderate control */
#define malloc m_alloc_d1
#define calloc c_alloc_d1
#define realloc re_alloc_d1
#define free fr_ee_d1
#endif

#ifdef DIAG2
/* and this to control all the memory constantly (slow) */
#define malloc m_alloc_d2
#define calloc c_alloc_d2
#define realloc re_alloc_d2
#define free fr_ee_d2
#endif

/*====================================================================*\
|| Prototypes der verschiedenen m_alloc-Funktionen, die bei der Suche ||
|| nach Fehlern in der Speicherverwaltung helfen koennen.             ||
|| Sie sind hier untergebracht, da sie vor dem einbinden von          ||
|| 'arith_const' bekannt sein muessen.                                ||
\*====================================================================*/
#ifdef __STDC__
void *m_alloc_d1(int );
void *c_alloc_d1(int , int );
void *re_alloc_d1(void *, int );
void fr_ee_d1(void *);

void *m_alloc_d2(int );
void *c_alloc_d2(int , int );
void *re_alloc_d2(void *, int );
void fr_ee_d2(void *);

void pointer_statistics(unsigned *,int);
#else
void *m_alloc_d1();
void *c_alloc_d1();
void *re_alloc_d1();
void fr_ee_d1();

void *m_alloc_d2();
void *c_alloc_d2();
void *re_alloc_d2();
void fr_ee_d2();

void pointer_statistics();
#endif

/*====================================================================*\
|| Pointer auf die verschiedenen m_alloc-Funktionen, die bei der Suche||
|| nach Fehlern in der Speicherverwaltung helfen koennen.             ||
\*====================================================================*/
/*
#ifdef __STDC__
extern int *(*m_alloc ) (int       );
extern int *(*c_alloc ) (int  , int);
extern int *(*re_alloc) (void *, int);
extern void (*fr_ee   ) (void *    );
#else
extern int *(*m_alloc ) (int       );
extern int *(*c_alloc ) (int  , int);
extern int *(*re_alloc) (void *, int);
extern void (*fr_ee   ) (void *    );
#endif
*/

#endif /* #ifndef _MALLOC_ */


