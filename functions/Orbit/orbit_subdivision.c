#include "typedef.h"
#include "utils.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: orbit_subdivision.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


typedef struct {
	int	dim;
	int	len;
	int	n;
	int	**v;
	} veclist;

static int *st_w;


static void 
vecmatmul (int *x, int **A, int *y, int dim)
{
	int	i, j, xi, *Ai;

	for (j = 0; j < dim; ++j)
		y[j] = 0;
	for (i = 0; i < dim; ++i)
	{
		if ((xi = x[i]) != 0)
		{
			Ai = A[i];
			for (j = 0; j < dim; ++j)
				y[j] += xi * Ai[j];
		}
	}
}

static int 
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

static int 
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
			low = search + 1;
		else if (cmp == -1)
			high = search - 1;
		else 
			return(sign * search);
	}
	return(sign * (V.n+low));
}



static void 
quicksort (	/*****	standard quicksort	*****/
    int **v,
    int inf,
    int sup,
    int dim
)
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

static void 
sortvecs (veclist *V)
{
	int	i, j;

	quicksort(V->v, 1, V->n, V->dim);
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


static int 
operate (int nr, int **A, veclist V)
{
	int	i, im;

	vecmatmul(V.v[abs(nr)], A, st_w, V.dim);
	if (nr < 0)
	{
		for (i = 0; i < V.dim; ++i)
			st_w[i] *= -1;
	}
/**********************************
        for(i=0;i<V.dim;i++)
          printf("%d  ", st_w[i]);
        printf("\n");
**********************************/
	im = numberof(st_w, V);
	if (abs(im) > V.n)
	{
		fprintf(stderr, "Error: image of vector %d not found\n", nr);
		exit (3);
	}
	return(im);
}




/**************************************************************************\
@---------------------------------------------------------------------------
@ int *orbit_subdivision(vecs, G, orbit_no)
@ matrix_TYP *vecs;
@ bravais_TYP *G;
@ int *orbit_no;
@
@ 'orbit_subdivision' calculates representatives of the orbit
@ of the group 'G' on the rows of the matrix 'vecs'
@ The action is assumed to be v -> vg^{tr} for v a row of 'vecs'
@ and g in 'G'.
@ It is assumed that -Identity is an element of 'G'
@ and for a row v the -v must not be contained as a row of 'vecs'.
@ Furthermore it is assumed, that the rows of 'vecs' are closed under
@ the action of 'G', so if vg^{tr} or -vg^{tr} is no row of 'vecs'
@ the programm exits (v rows of 'vecs'. g in 'G').
@ 
@ The number of orbits is returned via (int *orbit_no) and
@ the result is a pointer to integer of length 'orbit_no'
@ where the i-th entry of this pointer is the number of the
@ row of 'vecs' that is an representative of the i-th orbit.
@
@---------------------------------------------------------------------------
@
\**************************************************************************/
int *
orbit_subdivision (matrix_TYP *vecs, bravais_TYP *G, int *orbit_no)
{
	veclist	V;
	int	i, j, k, dim, nG, orblen, orbsum, orbnr, im, cnd;
	int	*orb, *flag, ***grp;
        int *Vvi;

	nG = G->gen_no;
        dim = G->gen[0]->cols;
	V.dim = dim;
        V.len = vecs->cols;
        V.n = vecs->rows;
	orb = (int*)xmalloc(V.n * sizeof(int));
	flag = (int*)xmalloc((V.n + 1) * sizeof(int));
	for (i = 0; i <= V.n; ++i)
		flag[i] = 0;
	st_w = (int*)xmalloc(V.dim * sizeof(int));
	V.v = (int**)xmalloc(((V.n)+1) * sizeof(int*));
        for(i=0;i<V.n;i++)
          V.v[i+1] = vecs->array.SZ[i];
        grp = (int ***)xmalloc(nG *sizeof(int **));
        for(i=0;i<nG;i++)
        {
          grp[i] = (int **)xmalloc(dim *sizeof(int *));
          for(j=0;j<dim;j++)
          {
            grp[i][j] = (int *)xmalloc(dim *sizeof(int));
            for(k=0;k<dim;k++)
              grp[i][j][k] = G->gen[i]->array.SZ[k][j];
          }
        }
	for (i = 1; i <= V.n; ++i)
	{
		Vvi = V.v[i];
		for (j = 0; j < V.dim  &&  Vvi[j] == 0; ++j);
		if (j < V.dim  &&  Vvi[j] < 0)
		{
			for (j = 0; j < V.dim; ++j)
				Vvi[j] *= -1;
		}
	}
	sortvecs(&V);
	for (i = 1; i <= V.n; ++i)
            vecs->array.SZ[i-1] = V.v[i];
	orbsum = 0;
	orbnr = 1;
	while (orbsum < 2*V.n)
	{
		for (i = 1; i <= V.n  &&  flag[i] != 0; ++i);
		orb[0] = i;
		flag[i] = orbnr;
		orblen = 1;
		cnd = 0;
		while (cnd < orblen)
		{
			for (i = 0; i < nG; ++i)
			{
				im = abs(operate(orb[cnd], grp[i], V));
				if (flag[im] == 0)
				{
					++orblen;
					orb[orblen-1] = im;
					flag[im] = orbnr;
				}
			}
			++cnd;
		}
/**********************************
		for (i = 0; i < V.dim; ++i)
			printf(" %2d", V.v[orb[0]][i]);
		printf("\t Length %d\n", 2*orblen);
**********************************/
		orbsum += 2*orblen;
		++orbnr;	
	}
        *orbit_no = orbnr-1;

        /* inserted the next 8 lines 15/1/97 */
        for (i=0;i<nG;i++){
           for (j=0;j<dim;j++){
              free(grp[i][j]);
           }
           free(grp[i]);
        }
        free(grp);
        free(V.v);

	free(orb);
        free(st_w);
	for (i = 1; i <= V.n; ++i)
           flag[i-1] = flag[i];
	return(flag);
}
