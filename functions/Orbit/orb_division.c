#include "typedef.h"
#include "matrix.h"
#include "sort.h"

/**************************************************************************\
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@ FILE: orb_division.c
@---------------------------------------------------------------------------
@---------------------------------------------------------------------------
@
\**************************************************************************/


/*****************************************
@ matrix_TYP *orbit_representatives(M, Manz, G, option, orbit_no, is_sorted)
@ matrix_TYP **M;
@ bravais_TYP *G;
@ int *option, Manz, *orbit_no;
@ int is_sorted;
@
@  The function 'orbit_representatives' calculates representatives
@  of a group 'G' (*bravais_TYP) on a set 'M' (**matrix_TYP) of matrices.
@  'Manz' denotes the number of matrices in 'M'.
@  The number of representatives is return via the pointer (int *orbit_no).
@  If the set 'M' is sorted (with respect to the order defined in the
@  function 'cmp_mat()') the function makes use of it by searching
@  orbit elements in a sorted list.
@  if is_sorted = 1, it is assumed that 'M' is sorted.
@                    CAUTION: This is not checked !
@                    So is_sorted = 1 for unsorted 'M' yields a wrong result.
@  if is_sorted = 0, it is assumed that 'M' is unsorted.
@   
@  The options are the same as in 'orbit_alg'.
******************************************/

matrix_TYP *orbit_representatives(M, Manz, G, option, orbit_no, is_sorted)
matrix_TYP **M;
bravais_TYP *G;
int *option, Manz, *orbit_no;
int is_sorted;
{
   matrix_TYP *erg, **or;
   int i,j,k, no;
   int *merk, found;

   extern matrix_TYP *init_mat();
   extern matrix_TYP **orbit_alg();
   extern int mat_comp();

  erg = init_mat(2, Manz, "");
  if((merk = (int *)calloc(Manz , sizeof(int))) == NULL)
  {
    printf("calloc of 'merk' in 'orbit_representatives' failed\n");
    exit(2);
  }

  i = 0;
  no = 0;
  while(i<Manz)
  {
    while(i<Manz && merk[i] != 0)
      i++;
    if(i<Manz)
    {
       erg->array.SZ[0][no] = i;
       or = orbit_alg(M[i], G, NULL, option, &(erg->array.SZ[1][no]));
       merk[i] = no+1;
       if(is_sorted == 1)
       {
         for(j=1;j<erg->array.SZ[1][no];j++)
         {
           found = mat_search(or[j], M, Manz, mat_comp);
           if(found != -1)
             merk[found] = no+1;
         }
       }
       else
       {
         for(j=1;j<erg->array.SZ[1][no];j++)
         {
           found = FALSE;
           for(k=i+1;k<Manz && found == FALSE;k++)
           {
             if(merk[k] != 1 && cmp_mat(or[j], M[k]) == 0)
             {
              merk[k] = no+1;
              found = TRUE;
             }
           }
         }
       }
       for(j=0;j<erg->array.SZ[1][no];j++)
         free_mat(or[j]);
       free(or);
       no++;
    }
  }
  real_mat(erg, 2, no);
  *orbit_no = no;
  free(merk);
  return(erg);
}
