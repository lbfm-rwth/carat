/*====================================================================*\
|| Prototypes der verschiedenen m_alloc-Funktionen, die bei der Suche ||
|| nach Fehlern in der Speicherverwaltung helfen koennen.             ||
|| Sie sind hier untergebracht, da sie vor dem einbinden von          ||
|| 'arith_const' bekannt sein muessen.                                ||
\*====================================================================*/
int *m_alloc_d1(int );
int *c_alloc_d1(int , int );
int *re_alloc_d1(int *, int );
void fr_ee_d1(int *);

int *m_alloc_d2(int );
int *c_alloc_d2(int , int );
int *re_alloc_d2(int *, int );
void fr_ee_d2(int *);
