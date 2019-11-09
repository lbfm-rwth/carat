#include "typedef.h"
#include "sort.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: quicksort.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/
/**********************************************************************\
@	standard quicksort algorithms
@   the are called like 'quicksort(V, 0, n-1, comp)
@   where V is a pointer to n objects and comp is a function that
@   compares the objects.
\**********************************************************************/




/**************************************************************************\
@---------------------------------------------------------------------------
@ void mat_quicksort(M, inf, sup, comp)
@ matrix_TYP **M;
@ int	inf, sup, (*comp)();
@
@ sorts a list of matrices M from M[inf] to M[sup] with respect to comp.
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
mat_quicksort (matrix_TYP **M, int inf, int sup, int (*comp)(const matrix_TYP *, const matrix_TYP *))
{
	int	low, med, high;
        matrix_TYP *tmp;

	if(inf+1 < sup)
	{
/* interchange v[inf] and v[med] to avoid worst case for pre-sorted lists */
		med = (inf+1+sup)/2;
		tmp = M[inf];
		M[inf] = M[med];
		M[med] = tmp;
		low = inf+1;
		high = sup;
		while(low < high)
		{
			while(comp(M[inf], M[low]) >= 0  &&  low < sup)
				++low;
			while(comp(M[inf], M[high]) <= 0  &&  high > inf)
				--high;
			if(low < high)
			{
				tmp = M[high];
				M[high] = M[low];
				M[low] = tmp;
			}
		}
		if(inf != high)
		{
			tmp = M[high];
			M[high] = M[inf];
			M[inf] = tmp;
		}
		mat_quicksort(M, inf, high-1, comp);
		mat_quicksort(M, high+1, sup, comp);
	}
	if(inf+1 == sup)
	{
		if(comp(M[inf], M[sup]) == 1)
		{
			tmp = M[inf];
			M[inf] = M[sup];
			M[sup] = tmp;
		}
	}
}

/**************************************************************************\
@---------------------------------------------------------------------------
@ void vec_quicksort(v, inf, sup, dim, comp)
@ int **v;
@ int	inf, sup, (*comp)(), dim;
@
@ sorts a list of vectors v from v[inf] to v[sup] with respect to comp.
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
vec_quicksort (int **v, int inf, int sup, int dim, int (*comp)(const int *, const int *, int))
{
	int	low, med, high;
        int *tmp;

	if(inf+1 < sup)
	{
/* interchange v[inf] and v[med] to avoid worst case for pre-sorted lists */
		med = (inf+1+sup)/2;
		tmp = v[inf];
		v[inf] = v[med];
		v[med] = tmp;
		low = inf+1;
		high = sup;
		while(low < high)
		{
			while(comp(v[inf], v[low], dim) >= 0  &&  low < sup)
				++low;
			while(comp(v[inf], v[high], dim) <= 0  &&  high > inf)
				--high;
			if(low < high)
			{
				tmp = v[high];
				v[high] = v[low];
				v[low] = tmp;
			}
		}
		if(inf != high)
		{
			tmp = v[high];
			v[high] = v[inf];
			v[inf] = tmp;
		}
		vec_quicksort(v, inf, high-1, dim, comp);
		vec_quicksort(v, high+1, sup, dim, comp);
	}
	if(inf+1 == sup)
	{
		if(comp(v[inf], v[sup], dim) == 1)
		{
			tmp = v[inf];
			v[inf] = v[sup];
			v[sup] = tmp;
		}
	}
}


/**************************************************************************\
@---------------------------------------------------------------------------
@ void pointer_mat_quicksort(v, inf, sup, rows, cols, comp)
@ int ***v;
@ int	inf, sup, rows, cols, (*comp)();
@
@---------------------------------------------------------------------------
@
@ sorts a list of 2-dimensional vectors v from v[inf] to v[sup]
@ with respect to comp.
\**************************************************************************/
void 
pointer_mat_quicksort (int ***v, int inf, int sup, int rows, int cols, int (*comp)(int **, int **, int, int))
{
	int	low, med, high;
        int **tmp;

	if(inf+1 < sup)
	{
/* interchange v[inf] and v[med] to avoid worst case for pre-sorted lists */
		med = (inf+1+sup)/2;
		tmp = v[inf];
		v[inf] = v[med];
		v[med] = tmp;
		low = inf+1;
		high = sup;
		while(low < high)
		{
			while(comp(v[inf], v[low], rows, cols) >= 0  &&  low < sup)
				++low;
			while(comp(v[inf], v[high], rows, cols) <= 0  &&  high > inf)
				--high;
			if(low < high)
			{
				tmp = v[high];
				v[high] = v[low];
				v[low] = tmp;
			}
		}
		if(inf != high)
		{
			tmp = v[high];
			v[high] = v[inf];
			v[inf] = tmp;
		}
		pointer_mat_quicksort(v, inf, high-1, rows, cols, comp);
		pointer_mat_quicksort(v, high+1, sup, rows, cols, comp);
	}
	if(inf+1 == sup)
	{
		if(comp(v[inf], v[sup], rows, cols) == 1)
		{
			tmp = v[inf];
			v[inf] = v[sup];
			v[sup] = tmp;
		}
	}
}
