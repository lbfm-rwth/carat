#ifdef __cplusplus
extern "C" {
#endif

/* check that this routine is included only once */
#ifndef _MALLOC_
#define _MALLOC_

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
void *m_alloc_d1(int );
void *c_alloc_d1(int , int );
void *re_alloc_d1(void *, int );
void fr_ee_d1(void *);

void *m_alloc_d2(int );
void *c_alloc_d2(int , int );
void *re_alloc_d2(void *, int );
void fr_ee_d2(void *);

void pointer_statistics(unsigned *,int);

#endif /* #ifndef _MALLOC_ */

#ifdef __cplusplus
}
#endif
