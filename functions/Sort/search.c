#include "typedef.h"
#include "sort.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: search.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
/**********************************************************************\
@	standard algorithms for searching vectors or matrices in
@       sortet  lists (sorted with respect to a function comp
@       that compares to elements in a list.
@
@       The functions return -1, if the element is not in the list
@       or the integer 0<= i < List_no where i denotes the
@       index where the element was found in the list.
\**********************************************************************/




/**************************************************************************\
@---------------------------------------------------------------------------
@ int mat_search(M, List, List_no, comp)
@ matrix_TYP *M, **List;
@ int	List_no, (*comp)();
@---------------------------------------------------------------------------
@
\**************************************************************************/
int 
mat_search (const matrix_TYP *M, matrix_TYP **List, int List_no, int (*comp)(const matrix_TYP *, const matrix_TYP *))
{
	int	low, med, high;
        int found = FALSE
           , test;

        low = 0; high = List_no-1;
        while(found == FALSE && low <= high)
        {
           med = (low + high)/2;
           test = comp(M, List[med]);
           if(test == 0)
             return(med);
           if(test > 0)
             low = med+1;
           else
             high = med-1;
        }
        return(-1);
}
