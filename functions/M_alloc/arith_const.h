/*====================================================================*\
|| Pointer auf die verschiedenen m_alloc-Funktionen, die bei der Suche||
|| nach Fehlern in der Speicherverwaltung helfen koennen.             ||
\*====================================================================*/

int *(*m_alloc ) (int        ) = (int *(*)(int        ))m_alloc_d1;
int *(*c_alloc ) (int   , int) = (int *(*)(int   , int))c_alloc_d1;
int *(*re_alloc) (void *, int) = (int *(*)(void *, int))re_alloc_d1;
void (*fr_ee   ) (void *     ) = (void (*)(void *     ))fr_ee_d1;

