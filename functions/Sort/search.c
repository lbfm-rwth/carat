#include"typedef.h"
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
int mat_search(M, List, List_no, comp)
matrix_TYP *M, **List;
int	List_no, (*comp)();
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



/**************************************************************************\
@---------------------------------------------------------------------------
@ int vec_search(M, List, List_no, dim, comp)
@ int *M, **List;
@ int	List_no, dim, (*comp)();
@---------------------------------------------------------------------------
@
\**************************************************************************/
int vec_search(M, List, List_no, dim, comp)
int *M, **List;
int	List_no, dim, (*comp)();
{
	int	low, med, high;
        int found = FALSE, test;

        low = 0; high = List_no-1;
        while(found == FALSE && low <= high)
        {
           med = (low + high)/2;
           test = comp(M, List[med], dim);
           if(test == 0)
             return(med);
           if(test == 1)
             low = med+1;
           else
             high = med-1;
        }
        return(-1);
}




/**************************************************************************\
@---------------------------------------------------------------------------
@ int pointer_mat_search(M, List, List_no, rows, cols, comp)
@ int **M, ***List;
@ int	List_no, rows, cols, (*comp)();
@---------------------------------------------------------------------------
@
\**************************************************************************/
int pointer_mat_search(M, List, List_no, rows, cols, comp)
int **M, ***List;
int	List_no, rows, cols, (*comp)();
{
	int	low, med, high;
        int found = FALSE, test;

        low = 0; high = List_no;
        while(found == FALSE && low <= high)
        {
           med = (low + high)/2;
           test = comp(M, List[med], rows, cols);
           if(test == 0)
             return(med);
           if(test == 1)
             low = med+1;
           else
             high = med-1;
        }
        return(-1);
}
