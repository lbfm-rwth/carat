/*****	This file includes some routines for 
	sorting lists and searching in sorted lists	*****/
#include "typedef.h"
#include "types.h"

/**********************************************************************\
|	compares the 1xn-vectors x and y lexicographically
|	returns 1 if x > y, 0 if x = y, -1 if x < y
\**********************************************************************/
int 
comp (int *x, int *y, int n)
{
	int	i;
  
	for (i = 0; i < n  &&  x[i] == y[i]; ++i);
	if (i == n)
		return(0);
	else if (x[i] > y[i])
		return(1);
	else
		return(-1);
}

/**********************************************************************\
|	searches for the vector vec in sorted 
|	list v := V.v between v[1] and v[V.n],
|	where the first non-zero entry in v[i] is > 0,
|	  returns i > 0 if vec = v[i],
|	  returns -i < 0 if -vec = v[i],
|	if vector is not found it returns
|	  V.n+i if v[i-1] < vec < v[i] or
|	  -(V.n+i) if v[i-1] < -vec < v[i]
|	if the return value is negative, vec is replaced by -vec
\**********************************************************************/
int 
numberof (int *vec, veclist V)
{
	int	i, sign, dim, low, high, search, cmp;

	sign = 1;
	dim = V.dim;
	low = 1;
	high = V.n;
	for (i = 0; i < dim  &&  vec[i] == 0; ++i);
	if (i < dim  &&  vec[i] < 0)
	{
		sign = -1;
		for (i = 0; i < dim; ++i)
			vec[i] *= -1;
	}
	while (low <= high)
	{
		search = low + (high-low)/2;
		cmp = comp(vec, V.v[search], dim);
		if (cmp == 1)
/* the vector is in the upper half */
			low = search + 1;
		else if (cmp == -1)
/* the vector is in the lower half */
			high = search - 1;
		else 
/* the vector is found */
			return(sign * search);
	}
/* if low > high the vector was not found */
	return(sign * (V.n+low));
}

/**********************************************************************\
|	sorts the V->n vectors v[1]...v[V->n] and 
|	deletes doublets, V->v is changed !!!
\**********************************************************************/
void 
sortvecs (veclist *V)
{
	int	i, j;

/* sort the vectors */
	quicksort(V->v, 1, V->n, V->dim);
/* now delete doublets */
	j = 1;
	for (i = 1; i+j <= V->n; ++i)
	{
		while (i+j <= V->n  &&  comp(V->v[i], V->v[i+j], V->dim) == 0)
		{
			free(V->v[i+j]);
			++j;
		}
		if (i+j <= V->n)
		{
			V->v[i+1] = V->v[i+j];
			if (comp(V->v[i], V->v[i+1], V->dim) == 1)
			{
				fprintf(stderr, "Error: v[%d] > v[%d]\n", i, i+1);
				exit (3);
			}
		}
	}
	V->n -= j-1;
}

/**********************************************************************\
|	standard quicksort
\**********************************************************************/
void 
quicksort (int **v, int inf, int sup, int dim)
{
	int	*tmp, low, med, high;

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
		quicksort(v, inf, high-1, dim);
		quicksort(v, high+1, sup, dim);
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
