#include"typedef.h"
#include "gmp.h"
/* #include "gmp-impl.h" */

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: MP_red_sort.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/**************************************************************************\
@---------------------------------------------------------------------------
@ void MP_reduction_sort(G,T,n)
@ MP_INT **G, **T;
@ int n;
@
@ The same as 'red_sort' in directory 'Reduction' for MP_INT** insted
@ of matrix_TYP*
@---------------------------------------------------------------------------
@
\**************************************************************************/
void 
MP_reduction_sort (MP_INT **G, MP_INT **T, int n)
{
  int i,j;
  int minpos;
  MP_INT min, *tmp, merk;

  mpz_init(&min);
  mpz_init(&merk);
 
  for(i=0;i<n;i++)
  {
    mpz_set(&min, &G[i][i]);
    minpos = i;
    for(j=i+1;j<n;j++)
    {
     if(mpz_cmp(&G[j][j], &min)<0)
     { mpz_set(&min, &G[j][j]); minpos = j;}
    }
    if(minpos != i)
    {
       tmp = G[i]; G[i] = G[minpos]; G[minpos] = tmp;
       tmp = T[i]; T[i] = T[minpos]; T[minpos] = tmp;
       for(j=0;j<n;j++)
       {
         mpz_set(&merk, &G[j][i]);
         mpz_set(&G[j][i], &G[j][minpos]);
         mpz_set(&G[j][minpos], &merk);
       }
    }
  }
  mpz_clear(&merk);
  mpz_clear(&min);
}
