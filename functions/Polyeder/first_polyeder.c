#include"typedef.h"
#include"polyeder.h"
#include"matrix.h"
#include"longtools.h"
/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: first_polyeder.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/




/**************************************************************************\
@---------------------------------------------------------------------------
@ polyeder_TYP *first_polyeder(mauern, anz)
@ wall_TYP **mauern;
@ int anz;
@
@ The arguments of 'first_polyeder' are linear inequalities
@ representent by 'mauern'.
@ 'anz' is the number of these equalities.
@ 'first_polyeder' returns NULL if there are not n = mauern[0]->dim
@ independent inequalities among the walls.
@ Otherwise, 'first_polyeder' selects n linear independent inequalities
@ and calculates the linear simplex defined by them.
@---------------------------------------------------------------------------
@
\**************************************************************************/
polyeder_TYP *
first_polyeder (wall_TYP **mauern, int anz)
{
  int i,j,k,n;
  int s, rrang, *take;
  matrix_TYP *M, *Mi;
  polyeder_TYP *erg;

  
  n = mauern[0]->dim;
  M = init_mat(anz, n, "");
  take = (int *)malloc(anz *sizeof(int));
  for(i=0;i<anz;i++)
      take[i] = 0;
  for(i=0;i<anz;i++)
    for(j=0; j<n; j++)
     M->array.SZ[i][j] = mauern[i]->gl[j];
/******************************
  put_mat(M, NULL, "matrix in first_polyeder", 0);
******************************/
  
  M->rows = 1;
  rrang = n;
  while(rrang != 0 && M->rows <= anz)
  {

    /* changed 19/8/97 tilman from:
    s = long_row_gauss(M); to: */
    s = long_row_basis(M,FALSE);

    s = n - s;
    if(s == rrang)
      take[M->rows-1] = FALSE;
    else
      take[M->rows-1] = TRUE;
    rrang = s;
    if(rrang != 0)
      M->rows++;
  }
  if(rrang != 0)
  {
    M->rows = anz;
    free_mat(M);
    free(take);
    return(NULL);
  }

  erg = init_polyeder(n, n);
  erg->is_closed = FALSE;
  erg->is_degenerate = FALSE;
  for(i=0; i<n;i++)
  {
    erg->vert[i] = init_vertex(n, n-1);
    erg->wall[i] = init_wall(n);
  }
  k = 0;
  M->rows = n;
  for(i=0; i<anz;i++)
  {
     if(take[i] == TRUE)
     {
        for(j=0;j<n;j++)
        {
          M->array.SZ[k][j] = mauern[i]->gl[j];
          erg->wall[k]->gl[j] = mauern[i]->gl[j];
        }
        erg->wall[k]->next_no = mauern[i]->next_no;  /* 7 lines anne 8/10/97 */
        erg->wall[k]->next = NULL;
        if(mauern[i]->next != NULL){		
           erg->wall[k]->next = (int**)malloc(mauern[i]->next_no*sizeof(int*));
           memcpy(erg->wall[k]->next, mauern[i]->next, 
                  mauern[i]->next_no * sizeof(int)); 
        }
        erg->wall[k]->mat = copy_mat(mauern[i]->mat); /*anne, 4.3.97(2 Zeile)*/
        erg->wall[k]->neu = mauern[i]->neu;
        /* inserted next 5 lines to correspond to first_fuber */
        erg->wall[k]->nproduct = mauern[i]->nproduct;
        erg->wall[k]->product = (int *) malloc(erg->wall[k]->nproduct *
                                               sizeof(int));
        for (j=0;j<erg->wall[k]->nproduct;j++)
               erg->wall[k]->product[j] = mauern[i]->product[j];

        k++;
     }
  }
/******************************
  put_mat(M, NULL, "matrix in first_polyeder", 0);
******************************/
  Mi = long_mat_inv(M);
  M->rows = anz;
  free_mat(M);
  free(take);
  for(i=0;i<n;i++)
  {
      k = 0;
      for(j=0;j<n;j++)
      {
         erg->vert[i]->v[j] = Mi->array.SZ[j][i];
         if(i!=j)
         { erg->vert[i]->wall[k] = j; k++;}
      }
  }
  free_mat(Mi);
  for(i=0;i<n;i++)
    normal_vertex(erg->vert[i]);
  return(erg);
}
