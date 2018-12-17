#include "typedef.h"
#include "matrix.h"
#include "polyeder.h"

fund_domain *first_fuber(mauern, anz)
wall_TYP **mauern;
int anz;
{
  int i,j,k,p;
  int s, n;

  fund_domain *erg;
  matrix_TYP *M, *A, *v, *tmp_v;
  int *take;
  int test;
  wall_TYP *vec;

  extern fund_domain *init_fund_domain();
  extern vertex_TYP *init_vertex_fuber();
  extern wall_TYP *init_wall_fuber();
  extern int row_gauss();
  extern matrix_TYP *solve_mat();
  extern matrix_TYP *init_mat();
  
  test = 0;
  n = mauern[0]->dim;
  M = init_mat(n, n, "");
  take = (int *)malloc(anz *sizeof(int));
  
  M->rows = 0;
  i=0;
  while(M->rows != n && i < anz)
  {
    for(j=0;j<n;j++)
      M->array.SZ[M->rows][j] = mauern[i]->gl[j];
    M->rows++;
    s = row_gauss(M);
    if(s == M->rows)
      take[i] = TRUE;
    else
    {
      take[i] = FALSE;
      M->rows--;
    }
    i++;
  }
  if(M->rows != n)
  {
    free_mat(M);
    free(take);
    return(NULL);
  }

  erg = init_fund_domain(n, n);
  erg->is_finite = FALSE;
  erg->is_closed = FALSE;
  for(i=0; i<n;i++)
  {
    erg->vert[i] = init_vertex_fuber(n, n-1);
    erg->wall[i] = init_wall_fuber(n);
  }
  k = 0;
  for(i=0;i<anz && k<n;i++)
  {
    if(take[i] == TRUE)
    {
     for(j=0;j<n;j++)
      erg->wall[k]->gl[j] = mauern[i]->gl[j];

     /* geaendert am 08.11.94jk: Information ueber die korrespondierende Matrix
        sollten uebernommen werden */
     erg->wall[k]->nproduct = mauern[i]->nproduct;
     erg->wall[k]->product = (int*) malloc (sizeof(int)*mauern[i]->nproduct);
     for (j=0;j<mauern[i]->nproduct;j++)
        erg->wall[k]->product[j] = mauern[i]->product[j];
     /* Ende der Aenderung */
     k++;
    }
  }
  for(i=0;i<n;i++)
  {
    k = 0;
    for(j=0;j<n;j++)
    {
      if(j != i)
      { erg->vert[i]->wall[k] = j; k++;}
    }
  }

  A = init_mat(n-1, n, "l");
  for(i=0;i<n;i++)
  {
    p=0;
    for(j=0; j<n;j++)
    {
      if(i != j)
      {
        for(k=0;k<n;k++)
          A->array.SZ[p][k] = erg->wall[j]->gl[k];
        p++;
      }
    }

    /* changed on 16/12/96 tilman. The function solve didn't exist
       anymore, so I used solve_mat */
    tmp_v = solve_mat(A);
    v = tr_pose(tmp_v);
    free_mat(tmp_v);

    for(j=0;j<n;j++)
      erg->vert[i]->v[j] = v->array.SZ[j][0];
    free_mat(v);
    test = wall_times_vertex_fuber(erg->wall[i], erg->vert[i]);
    if(test < 0)
    {
      for(j=0;j<n;j++)
        erg->vert[i]->v[j] = -erg->vert[i]->v[j];
    }
  }

  M->rows = n;
  free_mat(M);
  free_mat(A);
  free(take);
  return(erg);
}
/*{{{}}}*/
